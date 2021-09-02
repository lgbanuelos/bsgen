import copy
import math
import multiprocessing
import os
import pandas as pd
import random
import re
import subprocess


from functools import reduce

from datetime import datetime
from datetime import timedelta
from pathlib import Path

from time import process_time

from pm4py.algo.conformance.alignments.petri_net import algorithm as alignments
from pm4py.algo.discovery.inductive import algorithm as inductive_miner
from pm4py.algo.filtering.log.variants import variants_filter
from pm4py.algo.simulation.playout.petri_net import algorithm as simulator
from pm4py.objects.log.exporter.xes import exporter as xes_exporter
from pm4py.objects.log.obj import EventLog
from pm4py.objects.conversion.log import converter as log_converter
from pm4py.objects.petri_net.exporter.variants import pnml
from pm4py.objects.petri_net.importer import importer as pnml_importer

from pm4py.objects.petri_net.utils.final_marking import discover_final_marking
from pm4py.objects.petri_net.utils.initial_marking import discover_initial_marking

from pm4py.util import constants


def inductive(log):
    return inductive_miner.apply(log)

def compute_number_of_fitting_traces(model, log):
    net, initial_marking, final_marking = model
    # initial_marking = discover_initial_marking(net)
    # final_marking = discover_final_marking(net)
    # I am assuming that the net is handling properly the silent transitions
    aligned_traces = alignments.apply_log(log, net, initial_marking, final_marking)
    fitting = reduce(lambda acc, d: acc + 1 if math.isclose(d['fitness'], 1.0) else acc, aligned_traces, 0)
    return fitting

# def compute_goodness(model, bslog): #, ouniqtraces):
#     net, initial_marking, final_marking = model
#     (bsuniqlog, uniqtraces) = dedup(bslog)
#     # I am assuming that the net is handling properly the silent transitions
#     aligned_traces = alignments.apply_log(bsuniqlog, net, initial_marking, final_marking)
#     fitting = reduce(lambda acc, d: acc + 1 if math.isclose(d['fitness'], 1.0) else acc, aligned_traces, 0)

#     return len(bsuniqlog), fitting, fitting/len(bsuniqlog)

    # initial_marking = discover_initial_marking(net)
    # final_marking = discover_final_marking(net)

def simulate_petri_net(S, number_of_traces):
    net, _, _ = pnml_importer.apply(str(S))
    initial_marking = discover_initial_marking(net)
    print(initial_marking)
    
    # Fix silent transitions
    for t in net.transitions:
        if t.label == t.name:
            t.label = None
        
    return simulator.apply(net, initial_marking, variant=simulator.Variants.BASIC_PLAYOUT, parameters={simulator.Variants.BASIC_PLAYOUT.value.Parameters.NO_TRACES: number_of_traces})

def export_log(log, folder_name, system_name, qualifier, n):
    filename = os.path.join(folder_name, f'{system_name}_{qualifier}_{n}.xes')
    xes_exporter.apply(log, filename)
    return filename

def export_gzipped_log(log, folder_name, system_name, qualifier, n):
    filename = os.path.join(folder_name, f'{system_name}_{qualifier}_{n}.xes.gz')
    xes_exporter.apply(log, filename, parameters={"compress": True})
    return filename

def export_net(model, folder_name, system_name, qualifier, n):
    net, initial_marking, final_marking = model
    
    # Convert Petri net into an eTree (in memory PNML representation)
    tree = pnml.export_petri_tree(net, initial_marking, final_marking)
    for e in tree.findall(".//initialMarking/.."):
        e.remove(e.find("initialMarking"))
    # Remove name element for all "invisible" transitions
    for e in tree.findall(".//*[@activity = '$invisible$']/.."):
        e.remove(e.find('name'))
        
    # Store the resulting PNML in a file
    filename = os.path.join(folder_name, f'{system_name}_{qualifier}_{n}.pnml')    
    tree.write(filename, pretty_print=True, xml_declaration=True, encoding=constants.DEFAULT_ENCODING)
    return filename

def augment_events(log_as_list, i = 1):
    newlist = []
    currenttime = datetime.now()
    for t in log_as_list:
        newt = []
        for e in t:
            e['case'] = str(i)
            if 'orgcase' not in e:
                e['orgcase'] = str(i)
            e['time:timestamp'] = currenttime
            currenttime += timedelta(seconds=1)
            newt.append(e)
        newlist.append(newt)
        i += 1
    return newlist

def augment_events_prime(log_as_list, i = 1):
    newlist = []
    currenttime = datetime.now()
    for t in log_as_list:
        for e in t:
            e['case'] = str(i)
            if 'orgcase' not in e:
                e['orgcase'] = str(i)
            e['time:timestamp'] = currenttime
            currenttime += timedelta(seconds=1)
            newlist.append(e)
        i += 1
    return newlist


def find_site_set(str1, str2, n):
    sites = []
    for p1 in range(len(str1) - n + 1):
        for p2 in range(len(str2) - n + 1):
            if str1[p1:p1+n] == str2[p2:p2+n]:
                sites.append((p1,p2))
    return sites

def crossover(t1, t2, n):
    labels1 = [e['concept:name'] for e in t1]
    labels2 = [e['concept:name'] for e in t2]
    sites = find_site_set(labels1, labels2, n)
    if sites == []:
        return (0, 0)
    else:
        return random.choice(sites)

def log_estimation_with_breeding(L, Ln, l, p = 0.5):
    nL = []
    while len(nL) < len(L):
        t1 = random.choice(L)
        t2 = random.choice(Ln)
        if random.random() < p:
            (x, y) = crossover(t1, t2, l)
            if x < 1:
                nL.append(copy.deepcopy(t1))
                nL.append(copy.deepcopy(t2))
            else:
                t3 = copy.deepcopy(t1[:x + l] + t2[y + l:])
                t4 = copy.deepcopy(t2[:y + l] + t1[x + l:])
                nL.append(t3)
                nL.append(t4)
        else:
            nL.append(copy.deepcopy(t1))
            nL.append(copy.deepcopy(t2))
    return nL

def dedup(log):
    variants = variants_filter.get_variants(log)
    unique_traces = []
    for key in variants.keys():
        unique_traces.append(variants[key][0])
    new_log = EventLog(attributes=log.attributes, extensions=log.extensions, globals=log._omni,
                       classifiers=log.classifiers, properties=log.properties)
    new_log._list = unique_traces
    return new_log

def sample_with_replacement(log, no_traces):
    traces = [copy.deepcopy(random.choice(log)) for i in range(no_traces)]
    parameters = {log_converter.Variants.TO_EVENT_LOG.value.Parameters.CASE_ID_KEY: 'case'}
    df = pd.DataFrame(augment_events_prime(traces))
    new_log = log_converter.apply(df, parameters=parameters, variant=log_converter.Variants.TO_EVENT_LOG)
    return new_log

def log_sample_with_breeding(log, number_of_generations, number_of_traces, k, p):
    log_as_list = []
    for trace in log:
        log_as_list.append([copy.deepcopy(e)  for e in trace])
    logs = copy.deepcopy(log_as_list)
    logn = copy.deepcopy(log_as_list)
    overall_size = len(log)    
    for _ in range(number_of_generations):
        new_generation = log_estimation_with_breeding(log_as_list, logn, k, p)
        logn = augment_events(new_generation, overall_size + 1)
        logs += logn
        overall_size = len(logs)
        
    return sample_with_replacement(logs, number_of_traces)

def entropia_precision(modelFilename, sysOrLogFilename): #, MAX_TIME = 300):
    stepTwo = subprocess.check_output(['java', '-jar', 'jbpt-pm-entropia-1.6.jar', '-emp', '-t', f"-rel={sysOrLogFilename}",  f"-ret={modelFilename}"]) #, timeout=MAX_TIME)
    stepTwo = stepTwo.decode("utf-8").split('\n')
    print(stepTwo[-1])
    return float(re.search("\d+\.\d+", stepTwo[-1]).group())

def entropia_recall(modelFilename, sysOrLogFilename):
    stepOne = subprocess.check_output(['java', '-jar', 'jbpt-pm-entropia-1.6.jar', '-emr', '-t', f"-rel={sysOrLogFilename}",  f"-ret={modelFilename}"]) #, timeout=MAX_TIME)
    stepOne = stepOne.decode("utf-8").split('\n')
    print(stepOne[-1])
    return float(re.search("\d+\.\d+", stepOne[-1]).group())

def entropia_coverage(modelFilename, sysOrLogFilename): #, MAX_TIME = 300):
    stepOne = subprocess.check_output(['java', '-jar', 'jbpt-pm-entropia-1.6.jar', '-emr', '-t', f"-rel={sysOrLogFilename}",  f"-ret={modelFilename}"]) #, timeout=MAX_TIME)
    stepOne = stepOne.decode("utf-8").split('\n')
    print(stepOne[-1])
    stepTwo = subprocess.check_output(['java', '-jar', 'jbpt-pm-entropia-1.6.jar', '-emp', '-t', f"-rel={sysOrLogFilename}",  f"-ret={modelFilename}"]) #, timeout=MAX_TIME)
    stepTwo = stepTwo.decode("utf-8").split('\n')
    print(stepTwo[-1])
    # try:
    cSM = float(re.search("\d+\.\d+", stepOne[-1]).group())
    cMS = float(re.search("\d+\.\d+", stepTwo[-1]).group())
    # except:
    #     cSM = cMS = 1
    return cSM, cMS

def BootstrapGeneralizationDataCollection(S, LS, LSp, LSM, PDT, EN, LEM, GM, K, p, NoG, outdir = 'output'):
    # Set up
    system_name = S.stem
    folder_name = os.path.join(outdir, system_name)
    os.mkdir(folder_name)

    re1 = open(os.path.join(outdir, system_name, f"baseline.csv"), "a")
    re2 = open(os.path.join(outdir, system_name, f"bootstrap.csv"), "a")

    (precision, recall) = (1, 1)
    for (n, lsm) in [(n, lsm) for n in LS for lsm in LSM]:
        start_time = process_time()
        L = lsm(S, n)
        sampling_time = process_time() - start_time

        start_time = process_time()
        lname = export_log(L, folder_name, system_name, 'baseline', n)
        sample_log_exporting_time = process_time() - start_time

        for d in PDT:
            start_time = process_time()
            M = d(L)
            discovery_time = process_time() - start_time

            start_time = process_time()
            mname = export_net(M, folder_name, system_name, 'baseline', n)
            model_exporting_time = process_time() - start_time

            for g in GM:
                # try:
                start_time = process_time()
                ms_precision = entropia_precision(mname, str(S))
                ms_precision_time = process_time() - start_time

                start_time = process_time()
                ms_recall = entropia_recall(mname, str(S))
                ms_recall_time = process_time() - start_time
                # (ms_recall, ms_precision) = g(mname, str(S)) #, 300)


                # if ms_precision > 0.9: # or ms_recall > 0.9:
                #     continue
                start_time = process_time()
                ml_precision = entropia_precision(mname, lname)
                ml_precision_time = process_time() - start_time

                start_time = process_time()
                ml_recall = entropia_recall(mname, lname)
                ml_recall_time = process_time() - start_time
                # (ml_recall, ml_precision) = g(mname, lname) #, 300)

                re1.write(f"{system_name};{lname};{d.__name__};{g.__name__};{ms_precision};{ms_recall};{ml_precision};{ml_recall};{sample_with_replacement};{sample_log_exporting_time};{discovery_time};{model_exporting_time};{ms_precision_time};{ms_recall_time};{ml_precision_time};{ml_recall_time}\n")
                re1.flush()
                # except:
                #     print(f"TIMEOUT: Computation of MSP/MSR for '{system_name}' was stopped after 5 minutes!")

                for (np, m, lem, k, nog) in [(np, m, lem, k, nog) for np in LSp for m in EN for lem in LEM for k in K for nog in NoG]:
                    for i in range(m):
                        start_time = process_time()
                        Lstari = lem(L,nog,np,k,p)
                        bs_sampling_time = process_time() - start_time

                        start_time = process_time()
                        lstariname = export_log(Lstari, folder_name, system_name, f'star_{nog}_{k}_{p}_{i}', np)
                        bs_log_exporting_time = process_time() - start_time

                        start_time = process_time()
                        dlog = dedup(Lstari)
                        bs_log_dedup_time = process_time() - start_time

                        start_time = process_time()      
                        fitting = compute_number_of_fitting_traces(M, dlog)
                        fitness_time = process_time() - start_time

                        # export_gzipped_log(dlog, folder_name, system_name, f'dedup_{nog}_{k}_{p}_{i}', np)

                        start_time = process_time()
                        precision = entropia_precision(mname, lstariname)
                        precision_time = process_time() - start_time

                        start_time = process_time()
                        recall = entropia_recall(mname, lstariname)
                        recall_time = process_time() - start_time
                        # (recall, precision) = g(mname, lstariname)# , 1800)
                        re2.write(f"{system_name};{lstariname};{d.__name__};{g.__name__};")
                        re2.write(f"{k};{nog};{p};{np};{m};{i};{len(dlog)};{fitting};{precision};{recall};")
                        re2.write(f"{bs_sampling_time};{bs_log_exporting_time};{bs_log_dedup_time};{fitness_time};{precision_time};{recall_time}\n")
                        re2.flush()
                        # except:
                        #     print(f"TIMEOUT: Computation of MLP/MLR for '{system_name}' was stopped after 30 minutes!")
                        #     pass
                        os.remove(lstariname)
            re1.close()
            re2.close()
            return

if __name__ == "__main__":
    input_dir = '../data/icpm2020'
    output_dir = '../data/output'

    Systems =  [p for p in Path(input_dir).glob('**/PE*_04_*.pnml')]
    print(len(Systems))
    
    with multiprocessing.Pool() as pool:
        pool.starmap(BootstrapGeneralizationDataCollection, [(
            s,                          # One target system at a time 
            [100],                      # LS    Size of the log sample (Taken directly from the target system)
            [100,200,400,1000,2000,4000],              # LSp   Size of the bootstrapped sample
            [simulate_petri_net],       # LSM   Function to sample the target system (via simulation of the Petri net)
            [inductive],                # PDT   Algorithm for model discovery (Inductive miner)
            [100],                      # EN    Number of repetitions for bootstraping
            [log_sample_with_breeding], # LEM   Function to run trace breeding and sampling
            [entropia_coverage],        # GM    Function to compute both Model/System, Model/Log Precision & Recall 
            [2,3],                        # k     Size of the 
            1.0,                        # p     Breeding probability
            [100,200,400,1000,2000,4000],                      # NoG   Number of generations during breeding
            output_dir                  
            ) for s in Systems])

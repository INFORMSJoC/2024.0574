import numpy as np
from sim_anneal_paths import *
from compute_A_undirected import *
from compute_sols_undir import *
from compute_sols_undir_dw import *
from clean_cycles_sols import *
from graver_diff_sols import *
from edges_cost import *
import time
import argparse
import csv
from pathlib import Path
import os
import shutil

def load_data_graph(nodes_file, edges_file):

    read_nodes = []
    with open(nodes_file, 'r') as nodes_csv_file:
         reader = csv.reader(nodes_csv_file)
         for r in reader:
             read_nodes.append(r)
         nodes_csv_file.close()
    read_nodes = np.array(read_nodes[1:])

    data_edges = []
    with open(edges_file, 'r') as edges_csv_file:
         reader = csv.reader(edges_csv_file)
         for r in reader:
             data_edges.append(r)
         edges_csv_file.close()
    data_edges = np.int32(data_edges[1:])

    data_edges_dir = []
    cij_edges_dir = []
    tij_edges_dir = [] 
    lanes_edges_dir = []
    for i in range(data_edges.shape[0]):
        if (data_edges[i][0] < data_edges[i][1]):
           data_edges_dir.append([data_edges[i][0], data_edges[i][1]])
           cij_edges_dir.append(data_edges[i][3])
           tij_edges_dir.append(data_edges[i][2])
           lanes_edges_dir.append(data_edges[i][4])
    data_edges_dir = np.array(np.int32(data_edges_dir))
    cij_edges_dir = np.array(np.int32(cij_edges_dir))
    tij_edges_dir = np.array(np.int32(tij_edges_dir))
    lanes_edges_dir = np.array(np.int32(lanes_edges_dir))
    nodes_demand_size = np.int32(read_nodes[:, 2][np.int32(read_nodes[:, 4]) == 1])
    data_nodes = np.zeros((read_nodes.shape[0], 3))
    data_nodes[:,0] = np.arange(0, read_nodes.shape[0])
    data_nodes[:,1] = np.float64(read_nodes[:, 0])
    data_nodes[:,2] = np.float64(read_nodes[:, 1])
    nodes_demand = np.int32(data_nodes[:, 0][np.int32(read_nodes[:,4]) == 1])
    nodes_exit = np.int32(data_nodes[:, 0][np.int32(read_nodes[:,3]) == 1])
    return data_nodes, data_edges_dir, cij_edges_dir, tij_edges_dir, lanes_edges_dir, nodes_demand, nodes_exit, nodes_demand_size


def do_undir(cij_edges_dir, tij_edges_dir, lanes_edges_dir):

    cij_edges_undir = np.hstack([cij_edges_dir, cij_edges_dir])

    tij_edges_undir = np.hstack([tij_edges_dir, tij_edges_dir])

    lanes_edges_undir = np.hstack([lanes_edges_dir, lanes_edges_dir])

    return cij_edges_undir, tij_edges_undir, lanes_edges_undir


def do_feas_sols_cals(save_dir, data_nodes, data_edges_dir, nodes_demand, nodes_exit, n_samples):

    compute_sols_undir(save_dir, data_nodes, data_edges_dir, nodes_demand, nodes_exit, n_samples)


def get_sols_per_demand(save_dir, data_nodes, data_edges_dir, nodes_demand, nodes_exit, tij_edges_undir, A_mat_2n_noedges, near_exits=25):
    fsol_list_reduced = []
    for i_demand in range(nodes_demand.shape[0]):
        fi_l = np.genfromtxt(save_dir + '/feas_sols_sorted_' + str(i_demand) + '.txt', dtype=np.int32)
        node_demand_i = np.array([nodes_demand[i_demand]])
        b_2fr1, b_2fr1_noedges = compute_b_demand(data_nodes, node_demand_i, nodes_exit)
        #import pdb; pdb.set_trace()
        chk_i = (np.matmul(A_mat_2n_noedges, fi_l.T).T - b_2fr1_noedges.T == 0).all(axis=1)
        fi_l_i = fi_l[chk_i]
        fi_t_i = np.matmul(fi_l_i, tij_edges_undir)
        print(fi_l.shape, (~(chk_i)).sum(), fi_l_i.shape[0], fi_t_i.max(), fi_t_i[-1])
        if (len(fi_l_i) <= near_exits):
           fsol_list_reduced.append(fi_l_i)
        else:
           fsol_list_reduced.append(fi_l_i[0:near_exits])
    return fsol_list_reduced


def get_graver_basis(nodes_demand, fsol_list_reduced, graver_TIME):
    g_basis_reduced = []
    for i_demand in range(nodes_demand.shape[0]):
        t_ini = time.time()
        ker_sols = compute_kernel_feas_sols(fsol_list_reduced[i_demand])
        ker_sols_uniq = np.unique(ker_sols, axis=0)
        rand1 = np.arange(0, ker_sols_uniq.shape[0])
        np.random.shuffle(rand1)
        ker_sols_uniq_rand = ker_sols_uniq[rand1]
        g_basis_par_i = compute_gbasis_par(ker_sols_uniq_rand)
        g_basis_reduced.append(g_basis_par_i)
        t_fin = time.time()
        graver_TIME[i_demand] = t_fin - t_ini
        print(i_demand, len(g_basis_par_i))
    fsol_graver_r = np.vstack(g_basis_reduced)
    return g_basis_reduced, fsol_graver_r


def generate_random_feas(nodes_demand_size, tij_edges_undir, fsol_list_reduced):

    ini_sol_rand = np.zeros(tij_edges_undir.shape[0], dtype=np.int32)
    ini_fr_rand = np.zeros(tij_edges_undir.shape[0], dtype=np.int32)
    for i, fi_l in enumerate(fsol_list_reduced):
        for cc1 in range(nodes_demand_size[i]):
            ini_sol_rand += fi_l[np.random.randint(0, fi_l.shape[0])]
        ini_fr_rand += fi_l[np.random.randint(0, fi_l.shape[0])]
    return ini_sol_rand, ini_fr_rand


def generate_fixed_feas(nodes_demand_size, tij_edges_undir, fsol_list_reduced):
    ini_sol_rand = np.zeros(tij_edges_undir.shape[0], dtype=np.int32)
    ini_fr_rand = np.zeros(tij_edges_undir.shape[0], dtype=np.int32)
    for i, fi_l in enumerate(fsol_list_reduced):
        for cc1 in range(nodes_demand_size[i]):
            ini_sol_rand += fi_l[np.random.randint(0, fi_l.shape[0])]
        ini_fr_rand += fi_l[0]
    return ini_sol_rand, ini_fr_rand


def save_ini_sols(n_samples, nodes_demand_size, tij_edges_undir, fsol_list_reduced, save_dir, rand_ini_sols_dir='/rand_ini_sols_dir', rand_ini_fr_dir='/rand_ini_fr_dir'):
    for i in range(n_samples):
        ini_sol_rand, ini_fr_rand = generate_random_feas(nodes_demand_size, tij_edges_undir, fsol_list_reduced)
        ini_sol_rand.tofile(save_dir + '/' + rand_ini_sols_dir + '/ini_sol_rand_testfr' + str(i) + '.npy')
        ini_fr_rand.tofile(save_dir + '/' + rand_ini_fr_dir + '/ini_fr_rand_testfr' + str(i) + '.npy')
        print(i, n_samples)


def main_compute():

    parser = argparse.ArgumentParser()

    parser.add_argument('-nodes_file', '--nodes_file_cinput', type=str, help='input nodes file path')
    parser.add_argument('-edges_file', '--edges_file_cinput', type=str, help='input edges file path')
    parser.add_argument('-path_save_data', '--path_save_data_cinput', type=str, help='output directory path')
    parser.add_argument('-feas_sols', '--feas_sols_cinput', type=bool, default=False, help='compute Feasible solutions?')
    parser.add_argument('-n_samples', '--n_samples_cinput', type=int, default=10000, help='number of SA samples')
    parser.add_argument('-yens_k_shortest', '--yens_k_shortest_cinput', type=bool, default=False, help='run shortest path?')  
    parser.add_argument('-graver_basis', '--graver_basis_cinput', type=bool, default=False, help='compute Graver basis?')
    parser.add_argument('-dwave_sampling', '--dwave_sampling_cinput', type=bool, default=False, help='D-Wave sampling?')

    args = parser.parse_args()

    nodes_file                 = args.nodes_file_cinput
    edges_file                 = args.edges_file_cinput
    save_dir                   = args.path_save_data_cinput
    do_feas_sols_calc          = args.feas_sols_cinput
    n_samples                  = args.n_samples_cinput
    do_graver_basis_calc       = args.graver_basis_cinput
    do_yens_k_shortest         = args.yens_k_shortest_cinput
    do_DW                      = args.dwave_sampling_cinput

    n_demand = 100
    near_exits_k = 25
    
    if (not(Path(save_dir).is_dir())):
       os.makedirs(save_dir)
       os.mkdir(save_dir + '/rand_ini_sols_dir')
       os.mkdir(save_dir + '/rand_ini_fr_dir')
       os.mkdir(save_dir + '/sols_save_non0_optimal1')
       
    data_nodes, data_edges_dir, cij_edges_dir, tij_edges_dir, lanes_edges_dir, nodes_demand, nodes_exit, nodes_demand_size = load_data_graph(nodes_file, edges_file)

    cij_edges_dir[cij_edges_dir == 0] = 1
    np.savetxt(save_dir + '/data_nodes.txt', data_nodes)
    np.savetxt(save_dir + '/data_edges_dir.txt', data_edges_dir)
    np.savetxt(save_dir + '/cij_edges_dir.txt', cij_edges_dir)
    np.savetxt(save_dir + '/tij_edges_dir.txt', tij_edges_dir)
    np.savetxt(save_dir + '/lanes_edges_dir.txt', lanes_edges_dir)
    np.savetxt(save_dir + '/nodes_demand.txt', nodes_demand)
    np.savetxt(save_dir + '/nodes_exit.txt', nodes_exit)

    cij_edges_undir, tij_edges_undir, lanes_edges_undir = do_undir(cij_edges_dir, tij_edges_dir, lanes_edges_dir)
    np.savetxt(save_dir + '/cij_edges_undir.txt', cij_edges_undir)
    np.savetxt(save_dir + '/tij_edges_undir.txt', tij_edges_undir)
    np.savetxt(save_dir + '/lanes_edges_undir.txt', lanes_edges_undir)

    sols_SA_time = np.zeros(nodes_demand.shape[0])
    graver_TIME = np.zeros(nodes_demand.shape[0])
    t_ini = time.time()

    if do_feas_sols_calc:
       if do_yens_k_shortest:
          compute_sols_yens_k_shortest(save_dir, data_nodes, data_edges_dir, nodes_demand, nodes_exit, tij_edges_dir, n_samples, sols_SA_time)
       else:
          if do_DW:
              compute_sols_undir_outward_DW(save_dir, data_nodes, data_edges_dir, nodes_demand, nodes_exit, tij_edges_undir, n_samples, sols_SA_time)
          else:
              compute_sols_undir_outward(save_dir, data_nodes, data_edges_dir, nodes_demand, nodes_exit, tij_edges_undir, n_samples, sols_SA_time)
    A_mat_2n, A_mat_2n_noedges = compute_A_mat_undir(data_nodes, data_edges_dir, nodes_exit)
    fsol_list_reduced = get_sols_per_demand(save_dir, data_nodes, data_edges_dir, nodes_demand, nodes_exit, tij_edges_undir, A_mat_2n_noedges, near_exits=near_exits_k)
    for i_demand in range(len(fsol_list_reduced)):
        np.savetxt(save_dir + '/fsol_list_reduced_' + str(i_demand) + '.txt', fsol_list_reduced[i_demand])
    t_fin = time.time()
    time_SA = t_fin - t_ini
    
    t_ini = time.time()
    if do_graver_basis_calc:
       g_basis_reduced, fsol_graver_r = get_graver_basis(nodes_demand, fsol_list_reduced, graver_TIME)
       np.int32(fsol_graver_r).tofile(save_dir + '/fsol_graver_r.npy')
    else:
       fsol_graver_r = np.fromfile(save_dir + '/fsol_graver_r.npy', dtype=np.int32)
       fsol_graver_r = fsol_graver_r.reshape(-1, 2*data_edges_dir.shape[0])
    t_fin = time.time()
    time_graver = t_fin - t_ini
    np.savetxt(save_dir + '/times_sa_graver.txt', np.array([time_SA, time_graver]))
    ini_sol_rand, ini_fr_rand = generate_random_feas(nodes_demand_size, tij_edges_undir, fsol_list_reduced)
    np.int32(ini_sol_rand).tofile(save_dir + '/ini_sol_rand_testfr.npy')
    np.int32(ini_fr_rand).tofile(save_dir + '/ini_fr_rand_testfr.npy')
    n_rand_seeds = 1000
    save_ini_sols(n_rand_seeds, nodes_demand_size, tij_edges_undir, fsol_list_reduced, save_dir, rand_ini_sols_dir='/rand_ini_sols_dir', rand_ini_fr_dir='/rand_ini_fr_dir')
    
if __name__ == '__main__':
     main_compute()


## python main_gen_paths_random_instances.py --nodes_file_cinput='data/instances/random_instances/nodes_30_0.75_0.csv' --edges_file_cinput='data/instances/random_instances/edges_30_0.75_0.csv' --path_save_data_cinput='data/Q-HOPE_output/random_instances/graph_30_0.75_0' --nsamples_cinput=10000 --feas_sols_cinput=True --graver_basis_cinput=True

## python main_gen_paths_random_instances.py --nodes_file_cinput='data/instances/random_instances/nodes_10_0.75_0.csv' --edges_file_cinput='data/instances/random_instances/edges_10_0.75_0.csv' --path_save_data_cinput='data/Q-HOPE_output/random_instances/graph_10_0.75_0' --nsamples_cinput=1000 --feas_sols_cinput=True --graver_basis_cinput=True --dwave_sampling_cinput=True

## python main_gen_paths_random_instances.py --nodes_file_cinput='data/instances/random_instances/nodes_30_0.75_0.csv' --edges_file_cinput='data/instances/random_instances/edges_30_0.75_0.csv' --path_save_data_cinput='data/Q-HOPE_output/random_instances/graph_30_0.75_0' --feas_sols_cinput=True --graver_basis_cinput=True --yens_k_shortest_cinput=True


import numpy as np

dirname_from = './data/'
graph_dir_name0 = 'graph_'
graph_dir_name1 = 'graph_10000SA_'
graph_dir_name2 = 'graph_yens_k_shortest_'
graph_size = '20'
graph_den = '0.5'

print("Total time: SA+graver")
print(" ")
for i in range(12):
    try:
       a1_0 = np.fromfile(dirname_from + '/Q-HOPE_output/random_instances/' + graph_dir_name0 + graph_size + '_' + graph_den + '_' + str(i) + '/sols_save_non0_optimal1/e21sr_sol_srallFR_op1000.npy')
       b1_0 = np.fromfile(dirname_from + '/Q-HOPE_output/random_instances/' + graph_dir_name0 + graph_size + '_' + graph_den + '_' + str(i) + '/sols_save_non0_optimal1/e21time_sol_srallFR_op1000.npy')
       b1_0 = np.fromfile(dirname_from + '/Q-HOPE_output/random_instances/' + graph_dir_name0 + graph_size + '_' + graph_den + '_' + str(i) + '/sols_save_non0_optimal1/e21time_sol_srallFR_op1000.npy')
       a1_0l = np.fromfile(dirname_from + '/Q-HOPE_output/random_instances/' + graph_dir_name0 + graph_size + '_' + graph_den + '_' + str(i) + '/sols_save_non0_optimal1/e21sr_sol_nogama_op1000.npy')
       b1_0l = np.fromfile(dirname_from + '/Q-HOPE_output/random_instances/' + graph_dir_name0 + graph_size + '_' + graph_den + '_' + str(i) + '/sols_save_non0_optimal1/e21sr_time_nogama_op1000.npy')
       a1_0_chk = (a1_0 > 0)
       a1_0 = a1_0[a1_0_chk]
       b1_0 = b1_0[a1_0_chk]
       a1_0l = a1_0l[a1_0_chk]
       b1_0l = b1_0l[a1_0_chk]
       c1_0 = np.genfromtxt(dirname_from + '/Q-HOPE_output/random_instances/' + graph_dir_name0 + graph_size + '_' + graph_den + '_' + str(i) + '/times_sa_graver.txt')

       a1_1 = np.fromfile(dirname_from + '/Q-HOPE_output/random_instances/' + graph_dir_name1 + graph_size + '_' + graph_den + '_' + str(i) + '/sols_save_non0_optimal1/e21sr_sol_srallFR_op1000.npy')
       b1_1 = np.fromfile(dirname_from + '/Q-HOPE_output/random_instances/' + graph_dir_name1 + graph_size + '_' + graph_den + '_' + str(i) + '/sols_save_non0_optimal1/e21time_sol_srallFR_op1000.npy')
       a1_1l = np.fromfile(dirname_from + '/Q-HOPE_output/random_instances/' + graph_dir_name1 + graph_size + '_' + graph_den + '_' + str(i) + '/sols_save_non0_optimal1/e21sr_sol_nogama_op1000.npy')
       b1_1l = np.fromfile(dirname_from + '/Q-HOPE_output/random_instances/' + graph_dir_name1 + graph_size + '_' + graph_den + '_' + str(i) + '/sols_save_non0_optimal1/e21sr_time_nogama_op1000.npy')

       a1_1_chk = (a1_1 > 0)
       a1_1 = a1_1[a1_1_chk]
       b1_1 = b1_1[a1_1_chk]
       a1_1l = a1_1l[a1_1_chk]
       b1_1l = b1_1l[a1_1_chk]
       c1_1 = np.genfromtxt(dirname_from + '/Q-HOPE_output/random_instances/' + graph_dir_name1 + graph_size + '_' + graph_den + '_' + str(i) + '/times_sa_graver.txt')

       a1_2 = np.fromfile(dirname_from + '/Q-HOPE_output/random_instances/' + graph_dir_name2 + graph_size + '_' + graph_den + '_' + str(i) + '/sols_save_non0_optimal1/e21sr_sol_srallFR_op1000.npy')
       b1_2 = np.fromfile(dirname_from + '/Q-HOPE_output/random_instances/' + graph_dir_name2 + graph_size + '_' + graph_den + '_' + str(i) + '/sols_save_non0_optimal1/e21time_sol_srallFR_op1000.npy')
       b1_2l = np.fromfile(dirname_from + '/Q-HOPE_output/random_instances/' + graph_dir_name2 + graph_size + '_' + graph_den + '_' + str(i) + '/sols_save_non0_optimal1/e21sr_time_nogama_op1000.npy')
       a1_2l = np.fromfile(dirname_from + '/Q-HOPE_output/random_instances/' + graph_dir_name2 + graph_size + '_' + graph_den + '_' + str(i) + '/sols_save_non0_optimal1/e21sr_sol_nogama_op1000.npy')
       a1_2_chk = (a1_2 > 0)
       a1_2 = a1_2[a1_2_chk]
       b1_2 = b1_2[a1_2_chk]
       a1_2l = a1_2l[a1_2_chk]
       b1_2l = b1_2l[a1_2_chk]
       c1_2 = np.genfromtxt(dirname_from + '/Q-HOPE_output/random_instances/' + graph_dir_name2 + graph_size + '_' + graph_den + '_' + str(i) + '/times_sa_graver.txt')
       #b1_0l = np.array([0, 0])
       #b1_1l = np.array([0, 0])
       #b1_2l = np.array([0, 0])
       ##print(graph_size + " & " + graph_den + " & " + str(i) + " & ", "%.2fe+05" %(a1_0l.min()/1.0e5), " & ", "%.2f" %(b1_0.sum() + c1_0[0] + c1_0[1] + b1_0l.sum()), " & ","%.2fe+05" %(a1_1l.min()/1.0e5), " & ", "%.2f" %(b1_1.sum() + c1_1[0] + c1_1[1] + b1_1l.sum()), " & ", "%.2fe+05" %(a1_2l.min()/1.0e5), " & ", "%.2f" %(b1_2.sum() + c1_2[0] + c1_2[1] + b1_2l.sum()), " & & 14400", "\\")
       print(graph_size + " & " + graph_den + " & " + str(i) + " & ", "%.2e" %(a1_0l.min()), " & ", "%.2f" %(b1_0.sum() + c1_0[0] + c1_0[1] + b1_0l.sum()), " & ","%.2e" %(a1_1l.min()), " & ", "%.2f" %(b1_1.sum() + c1_1[0] + c1_1[1] + b1_1l.sum()), " & ", "%.2e" %(a1_2l.min()), " & ", "%.2f" %(b1_2.sum() + c1_2[0] + c1_2[1] + b1_2l.sum()), " & & 14400", "\\")
       #print(a1_1.argmin(), a1_1.argmin(), a1_1.argmin()) 
       #print(a1_0l.shape, a1_1l.shape, a1_2l.shape)
    except:
       pass
print(" ")



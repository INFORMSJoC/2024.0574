import numpy as np

def get_data(dirname, sols_type, id_instance, sols_dir, input_dir, do_sr1=True, pre_id=1):
    dirname1 = dirname + sols_type + id_instance + sols_dir
    a1_01 = np.fromfile(dirname1 + '/e21sr_sol_srallFR_op1000.npy')
    b1_0 = np.fromfile(dirname1 + '/e21time_sol_srallFR_op1000.npy')
    a1_0l = np.fromfile(dirname1 + '/e21sr_sol_nogama_op1000.npy')[0:]
    b1_0l = np.fromfile(dirname1 + '/e21sr_time_nogama_op1000.npy')[0:]
    a1_0 = a1_01[a1_01 > 0]
    b1_0 = b1_0[a1_01 > 0]
    a1_0l = a1_0l[a1_01 > 0]
    b1_0l = b1_0l[a1_01 > 0]
    c1_0 = np.genfromtxt(dirname + sols_type + id_instance + '/times_sa_graver.txt')

    return a1_0, b1_0, c1_0, a1_0l, b1_0l


def plot_normal(a1_0l, b1_0, c1_0, min0, label_line, id_sol_type, fig_id, b1_0l=None):

    t0i = c1_0[0] + c1_0[1]

    t0 = np.zeros(a1_0l.shape[0])
    min0obj = np.zeros(a1_0l.shape[0])
    for i0 in range(a1_0l.shape[0]):
        t0[i0] = t0i + b1_0[:i0+1].sum()
        if b1_0l is not None:
           t0[i0] = t0i + b1_0[:i0+1].sum() + b1_0l[:i0+1].sum()
        min0obj[i0] = a1_0l[:i0+1].min()

    matplotlib.pyplot.figure(fig_id)
    matplotlib.pyplot.plot(t0, min0obj/min0, label=label_line)
    matplotlib.pyplot.xlabel(r'$\mathrm{Time \ to \ solution}$ $(\mathrm{s})$', fontsize=20)
    matplotlib.pyplot.ylabel(r'$T_{evacuation}$/$T_{min}$', fontsize=20)
    matplotlib.pyplot.legend(loc=0)
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig('./fig_Tgraph_soltype_' + str(id_sol_type) + '.png')
    matplotlib.pyplot.show()


def plot_sorted(a1_0l, b1_0, c1_0, min0, label_line, id_sol_type, fig_id, b1_0l=None):

    t0i = c1_0[0] + c1_0[1]

    arg1 = np.argsort(a1_0l)[::-1]
    a1_0l = a1_0l[arg1]
    b1_0 = b1_0[arg1]
    b1_0l = b1_0l[arg1]

    t0 = np.zeros(a1_0l.shape[0])
    min0obj = np.zeros(a1_0l.shape[0])
    for i0 in range(a1_0l.shape[0]):
        t0[i0] = t0i + b1_0[:i0+1].sum()
        if b1_0l is not None:
           t0[i0] = t0i + b1_0[:i0+1].sum() + b1_0l[:i0+1].sum()
        min0obj[i0] = a1_0l[:i0+1].min()

    matplotlib.pyplot.figure(fig_id)
    matplotlib.pyplot.plot(t0, min0obj/min0, label=label_line)
    matplotlib.pyplot.xlabel(r'$\mathrm{Time \ to \ solution}$ $(\mathrm{s})$', fontsize=20)
    matplotlib.pyplot.ylabel(r'$T_{evacuation}$/$T_{min}$', fontsize=20)
    matplotlib.pyplot.legend(loc=0)
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig('./fig_Tgraph_soltype_' + str(id_sol_type) + 'w.png')
    matplotlib.pyplot.show()
    
dirname_Tgraph = './data/Q-HOPE_output/Tgraph'
dirname = dirname_Tgraph
input_dir = './data/instances/Tgraph_noNorm/'

t0_bnb = np.genfromtxt(dirname_Tgraph + '/bnb_logs/nodes_1_times.txt')
obj0_bnb = np.genfromtxt(dirname_Tgraph + '/bnb_logs/nodes_1_objective.txt')
t1_bnb = np.genfromtxt(dirname_Tgraph + '/bnb_logs/nodes_2_times.txt')
obj1_bnb = np.genfromtxt(dirname_Tgraph + '/bnb_logs/nodes_2_objective.txt')
t2_bnb = np.genfromtxt(dirname_Tgraph + '/bnb_logs/nodes_3_times.txt')
obj2_bnb = np.genfromtxt(dirname_Tgraph + '/bnb_logs/nodes_3_objective.txt')

min0 = 1.43e+06
min1 = 2.46e+07
min2 = 1.01e+08

matplotlib.pyplot.figure(0)
matplotlib.pyplot.plot(t0_bnb, obj0_bnb/min0, label='Branch & Bound')

matplotlib.pyplot.figure(1)
matplotlib.pyplot.plot(t1_bnb, obj1_bnb/min1, label='Branch & Bound')

matplotlib.pyplot.figure(2)
matplotlib.pyplot.plot(t2_bnb, obj2_bnb/min2, label='Branch & Bound')

matplotlib.pyplot.figure(3)
matplotlib.pyplot.plot(t0_bnb, obj0_bnb/min0, label='Branch & Bound')

matplotlib.pyplot.figure(4)
matplotlib.pyplot.plot(t1_bnb, obj1_bnb/min1, label='Branch & Bound')

matplotlib.pyplot.figure(5)
matplotlib.pyplot.plot(t2_bnb, obj2_bnb/min2, label='Branch & Bound')


sols_type = '/graph_10000SA1__'
sols_dir = '/Graver_walk'

a1_0, b1_0, c1_0, a1_0l, b1_0l = get_data(dirname, sols_type, '101', sols_dir, input_dir, do_sr1=False, pre_id=0)
a1_1, b1_1, c1_1, a1_1l, b1_1l = get_data(dirname, sols_type, '102', sols_dir, input_dir, do_sr1=False, pre_id=0)
a1_2, b1_2, c1_2, a1_2l, b1_2l = get_data(dirname, sols_type, '103', sols_dir, input_dir, do_sr1=False, pre_id=0)


print ("Normalized")
print('179' + " & " + '234' + " & ", "%.2fe+06" %(a1_0l[a1_0l>0].min()/1.0e6), " & ", "%.2f" %(b1_0.sum() + b1_0l.sum() + c1_0[0] + c1_0[1]), "\\")
print('179' + " & " + '234' + " & ", "%.2fe+07" %(a1_1l[a1_1l>0].min()/1.0e7), " & ", "%.2f" %(b1_1.sum() + b1_1l.sum() + c1_1[0] + c1_1[1]), "\\")
print('179' + " & " + '234' + " & ", "%.2fe+08" %(a1_2l[a1_2l>0].min()/1.0e8), " & ", "%.2f" %(b1_2.sum() + b1_2l.sum() + c1_2[0] + c1_2[1]), "\\")
print(a1_0l[a1_0l>0].shape, a1_1l[a1_1l>0].shape, a1_2l[a1_2l>0].shape)


plot_sorted(a1_0l, b1_0, c1_0, min0, 'Normalized', 0, 0, b1_0l=b1_0l)
plot_sorted(a1_1l, b1_1, c1_1, min1, 'Normalized', 1, 1, b1_0l=b1_1l)
plot_sorted(a1_2l, b1_2, c1_2, min2, 'Normalized', 2, 2, b1_0l=b1_2l)
plot_normal(a1_0l, b1_0, c1_0, min0, 'Normalized', 0, 3, b1_0l=b1_0l)
plot_normal(a1_1l, b1_1, c1_1, min1, 'Normalized', 1, 4, b1_0l=b1_1l)
plot_normal(a1_2l, b1_2, c1_2, min2, 'Normalized', 2, 5, b1_0l=b1_2l)

sols_dir = '/Graver_walk_tol'
a1_0, b1_0, c1_0, a1_0l, b1_0l = get_data(dirname, sols_type, '101', sols_dir, input_dir, do_sr1=False)
a1_1, b1_1, c1_1, a1_1l, b1_1l = get_data(dirname, sols_type, '102', sols_dir, input_dir, do_sr1=False)
a1_2, b1_2, c1_2, a1_2l, b1_2l = get_data(dirname, sols_type, '103', sols_dir, input_dir, do_sr1=False)

"""
assert np.isclose(a1_0l, t1_0l).all()
assert np.isclose(a1_1l, t1_1l).all()
assert np.isclose(a1_2l, t1_2l).all()
"""
print ("Normalized: Tol=1e-3")
print('179' + " & " + '234' + " & ", "%.2fe+06" %(a1_0l[a1_0l>0].min()/1.0e6), " & ", "%.2f" %(b1_0.sum() + b1_0l.sum() + c1_0[0] + c1_0[1]), "\\")
print('179' + " & " + '234' + " & ", "%.2fe+07" %(a1_1l[a1_1l>0].min()/1.0e7), " & ", "%.2f" %(b1_1.sum() + b1_1l.sum() + c1_1[0] + c1_1[1]), "\\")
print('179' + " & " + '234' + " & ", "%.2fe+08" %(a1_2l[a1_2l>0].min()/1.0e8), " & ", "%.2f" %(b1_2.sum() + b1_2l.sum() + c1_2[0] + c1_2[1]), "\\")
print(a1_0l[a1_0l>0].shape, a1_1l[a1_1l>0].shape, a1_2l[a1_2l>0].shape)

plot_sorted(a1_0l, b1_0, c1_0, min0, 'Normalized, Tol=1e-3', 0, 0, b1_0l=b1_0l)
plot_sorted(a1_1l, b1_1, c1_1, min1, 'Normalized, Tol=1e-3', 1, 1, b1_0l=b1_1l)
plot_sorted(a1_2l, b1_2, c1_2, min2, 'Normalized, Tol=1e-3', 2, 2, b1_0l=b1_2l)
plot_normal(a1_0l, b1_0, c1_0, min0, 'Normalized, Tol=1e-3', 0, 3, b1_0l=b1_0l)
plot_normal(a1_1l, b1_1, c1_1, min1, 'Normalized, Tol=1e-3', 1, 4, b1_0l=b1_1l)
plot_normal(a1_2l, b1_2, c1_2, min2, 'Normalized, Tol=1e-3', 2, 5, b1_0l=b1_2l)

sols_dir = '/Graver_walk_tol'
a1_0, b1_0, c1_0, a1_0l, b1_0l = get_data(dirname, sols_type, '1', sols_dir, input_dir)
#print_fr_unique(dirname + sols_type + '1' + sols_dir, 1, 100+1)
a1_1, b1_1, c1_1, a1_1l, b1_1l = get_data(dirname, sols_type, '2', sols_dir, input_dir)
#print_fr_unique(dirname + sols_type + '2' + sols_dir, 1, 100+1)
a1_2, b1_2, c1_2, a1_2l, b1_2l = get_data(dirname, sols_type, '3', sols_dir, input_dir)
#print_fr_unique(dirname + sols_type + '3' + sols_dir, 1, 100+1)
"""
assert np.isclose(a1_0, t1_0).all()
assert np.isclose(a1_1, t1_1).all()
assert np.isclose(a1_2, t1_2).all()
assert np.isclose(a1_0l, t1_0l).all()
assert np.isclose(a1_1l, t1_1l).all()
assert np.isclose(a1_2l, t1_2l).all()
"""
print("Unnormalized, Tol=1e-3")
print('179' + " & " + '234' + " & ", "%.2fe+06" %(a1_0l[a1_0l>0].min()/1.0e6), " & ", "%.2f" %(b1_0.sum() + b1_0l.sum() + c1_0[0] + c1_0[1]), "\\")
print('179' + " & " + '234' + " & ", "%.2fe+07" %(a1_1l[a1_1l>0].min()/1.0e7), " & ", "%.2f" %(b1_1.sum() + b1_1l.sum() + c1_1[0] + c1_1[1]), "\\")
print('179' + " & " + '234' + " & ", "%.2fe+08" %(a1_2l[a1_2l>0].min()/1.0e8), " & ", "%.2f" %(b1_2.sum() + b1_2l.sum() + c1_2[0] + c1_2[1]), "\\")
print(a1_0l[a1_0l>0].shape, a1_1l[a1_1l>0].shape, a1_2l[a1_2l>0].shape)

plot_sorted(a1_0l, b1_0, c1_0, min0, 'Unnormalized, Tol=1e-3', 0, 0, b1_0l=b1_0l)
plot_sorted(a1_1l, b1_1, c1_1, min1, 'Unnormalized, Tol=1e-3', 1, 1, b1_0l=b1_1l)
plot_sorted(a1_2l, b1_2, c1_2, min2, 'Unnormalized, Tol=1e-3', 2, 2, b1_0l=b1_2l)
plot_normal(a1_0l, b1_0, c1_0, min0, 'Unnormalized, Tol=1e-3', 0, 3, b1_0l=b1_0l)
plot_normal(a1_1l, b1_1, c1_1, min1, 'Unnormalized, Tol=1e-3', 1, 4, b1_0l=b1_1l)
plot_normal(a1_2l, b1_2, c1_2, min2, 'Unnormalized, Tol=1e-3', 2, 5, b1_0l=b1_2l)

nseed0 = 100
nseed1 = 85
nseed2 = 43

sols_dir = '/Graver_walk'
nseed0 = 100 
a1_0, b1_0, c1_0, a1_0l, b1_0l = get_data(dirname, sols_type, '1', sols_dir, input_dir)
a1_0i, b1_0i, c1_0i, a1_0il, b1_0il = a1_0[0:nseed0], b1_0[0:nseed0], c1_0[0:nseed0], a1_0l[0:nseed0], b1_0l[0:nseed0]

a1_1, b1_1, c1_1, a1_1l, b1_1l = get_data(dirname, sols_type, '2', sols_dir, input_dir)
a1_1i, b1_1i, c1_1i, a1_1il, b1_1il = a1_1[0:nseed1], b1_1[0:nseed1], c1_1[0:nseed1], a1_1l[0:nseed1], b1_1l[0:nseed1]

a1_2, b1_2, c1_2, a1_2l, b1_2l = get_data(dirname, sols_type, '3', sols_dir, input_dir)
a1_2i, b1_2i, c1_2i, a1_2il, b1_2il = a1_2[0:nseed2], b1_2[0:nseed2], c1_2[0:nseed2], a1_2l[0:nseed2], b1_2l[0:nseed2]

"""
assert np.isclose(a1_0i, t1_0i).all()
assert np.isclose(a1_1i, t1_1i).all()
assert np.isclose(a1_2i, t1_2i).all()
assert np.isclose(a1_0il, t1_0il).all()
assert np.isclose(a1_1il, t1_1il).all()
assert np.isclose(a1_2il, t1_2il).all()
"""
print("Unnormalized")
print('179' + " & " + '234' + " & ", "%.2fe+06" %(a1_0il.min()/1.0e6), " & ", "%.2f" %(b1_0i.sum() + c1_0i[0] + c1_0i[1] + b1_0il.sum()), "\\")
print('179' + " & " + '234' + " & ", "%.2fe+07" %(a1_1il.min()/1.0e7), " & ", "%.2f" %(b1_1i.sum() + c1_1i[0] + c1_1i[1] + b1_1il.sum()), "\\")
print('179' + " & " + '234' + " & ", "%.2fe+08" %(a1_2il.min()/1.0e8), " & ", "%.2f" %(b1_2i.sum() + c1_2i[0] + c1_2i[1] + b1_2il.sum()), "\\")
print(a1_0i.shape, a1_1i.shape, a1_2i.shape)

matplotlib.pyplot.figure(0)
matplotlib.pyplot.scatter(b1_0i.sum() + c1_0i[0] + c1_0i[1] + b1_0il.sum(), a1_0il.min()/min0, marker='*', color='purple')
plot_sorted(a1_0il, b1_0i, c1_0i, min0, r'Unnormalized ($^{*}$M=' + str(nseed0) + ')', 0, 0, b1_0l=b1_0il)
matplotlib.pyplot.figure(1)
matplotlib.pyplot.scatter(b1_1i.sum() + c1_1i[0] + c1_1i[1] + b1_1il.sum(), a1_1il.min()/min1, marker='*', color='purple')
plot_sorted(a1_1il, b1_1i, c1_1i, min1, r'Unnormalized ($^{*}$M=' + str(nseed1) + ')', 1, 1, b1_0l=b1_1il)
matplotlib.pyplot.figure(2)
matplotlib.pyplot.scatter(b1_2i.sum() + c1_2i[0] + c1_2i[1] + b1_2il.sum(), a1_2il.min()/min2, marker='*', color='purple')
plot_sorted(a1_2il, b1_2i, c1_2i, min2, r'Unnormalized ($^{*}$M=' + str(nseed2) + ')', 2, 2, b1_0l=b1_2il)
matplotlib.pyplot.figure(3)
matplotlib.pyplot.scatter(b1_0i.sum() + c1_0i[0] + c1_0i[1] + b1_0il.sum(), a1_0il.min()/min0, marker='*', color='purple')
plot_normal(a1_0il, b1_0i, c1_0i, min0, r'Unnormalized ($^{*}$M=' + str(nseed0) + ')', 0, 3, b1_0l=b1_0il)
matplotlib.pyplot.figure(4)
matplotlib.pyplot.scatter(b1_1i.sum() + c1_1i[0] + c1_1i[1] + b1_1il.sum(), a1_1il.min()/min1, marker='*', color='purple')
plot_normal(a1_1il, b1_1i, c1_1i, min1, r'Unnormalized ($^{*}$M=' + str(nseed1) + ')', 1, 4, b1_0l=b1_1il)
matplotlib.pyplot.figure(5)
matplotlib.pyplot.scatter(b1_2i.sum() + c1_2i[0] + c1_2i[1] + b1_2il.sum(), a1_2il.min()/min2, marker='*', color='purple')
plot_normal(a1_2il, b1_2i, c1_2i, min2, r'Unnormalized ($^{*}$M=' + str(nseed2) + ')', 2, 5, b1_0l=b1_2il)

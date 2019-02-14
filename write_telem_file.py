import sys

def writting(X, tau, t_max):
    original = sys.stdout
    ouf = open('D:\\redmine_svn\\497\\work\\progs\\telemetr_errors_bins.txt', 'w')
    sys.stdout = ouf
    print 'alpha', 'betta', 'ksi', 'inc_Vxg', 'inc_Vyg', 'inc_Vzg', 'inc_phi', 'inc_lambda', 'inc_h'
    for t in range(int(t_max/tau)):
        print X[t][0], X[t][1], X[t][2], X[t][3], X[t][4], X[t][5], X[t][6], X[t][7], X[t][8]
       #ouf.writelines(line)
    sys.stdout
    ouf.close()

def writting2(X, tau, t_max):
    ouf = open('D:\\redmine_svn\\497\\work\\progs\\telemetr_errors_bins.txt', 'w')
    shet = 0
    ouf.write('T_mod ')
    ouf.write('alpha ')
    ouf.write('betta ')
    ouf.write('ksi ')
    ouf.write('d_Vxg ')
    ouf.write('d_Vyg ')
    ouf.write('d_Vzg ')
    ouf.write('d_phi ')
    ouf.write('d_lambda ')
    ouf.write('d_h ')
    ouf.write('d_S\n')
    for t in range(int(t_max/tau)):
        shet = shet + tau
        ouf.write(str(shet))
        ouf.write(' ')
        ouf.write(str(X[t][0]*57.3))
        ouf.write(' ')
        ouf.write(str(X[t][1]*57.3))
        ouf.write(' ')
        ouf.write(str(X[t][2]*57.3))
        ouf.write(' ')
        ouf.write(str(X[t][3]))
        ouf.write(' ')
        ouf.write(str(X[t][4]))
        ouf.write(' ')
        ouf.write(str(X[t][5]))
        ouf.write(' ')
        ouf.write(str(X[t][6]))
        ouf.write(' ')
        ouf.write(str(X[t][7]))
        ouf.write(' ')
        ouf.write(str(X[t][8]))
        ouf.write(' ')
        ouf.write(str((X[t][6]*9.8)/pow(1.24*1e-3, 2)))
        ouf.write('\n')
    ouf.close()
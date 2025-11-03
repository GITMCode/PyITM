from pyitm.general import system

system.run_command('bin/plot_logfile.py data/log00000002.dat data/log00000154.dat -vars 0 17 18 19 20 21 22 -plotfile twologfilesplotted')


system.run_command('bin/plot_logfile.py data/log00000002.dat -vars 0 17 18 19 20 21 22 -plotfile onelogfile')

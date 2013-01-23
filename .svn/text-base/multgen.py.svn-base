import os

generations = [50,50,50,50]
name = 'millremove'
runs = 1

for gen in generations:
    filename = name + '%02d_%d.txt' % (runs, gen)
    plotfile = name + '%02d_%d_status.txt' % (runs, gen)
    command = 'python genalg.py ' + str(gen) + ' ' + filename + ' ' + plotfile
    os.system(command)
    runs += 1

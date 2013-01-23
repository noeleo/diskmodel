import matplotlib.pyplot as plt
import sys

def process(fileName):
    # read in data from file
    history = []
    f = open(fileName)
    i = 0
    curGen = []
    for line in f:
        data = map(float, line.split())
        curGen.append(data)
        i += 1
        if i == 100:
            history.append(curGen)
            curGen = []
            i = 0
    
    # find averages
    averages = [[] for x in range(9)]
    for generation in history:
        local_sum = [0]*9
        # find local sum
        for model in generation:
            for param in range(9):
                local_sum[param] += model[param]
        # divide by number of members in population to get average
        for param in range(9):
            averages[param].append(local_sum[param]/100.0)
            
    # plot the average chi values over generations
    plt.figure(1)
    # chis
    """
    plt.subplot(221)
    plt.plot(averages[8])
    plt.xlabel('Generation')
    plt.ylabel('Average Chi')
    plt.title('Chi-Squared (not reduced)')
    plt.subplot(223)
    plt.plot(averages[6])
    plt.xlabel('Generation')
    plt.ylabel('Average Chi')
    plt.title('SED')
    plt.subplot(224)
    plt.plot(averages[7])
    plt.xlabel('Generation')
    plt.ylabel('Average Chi')
    plt.title('Visibilities')
    """
    # DOES NOT WORK FOR NON-REDUCED
    #plt.ylim(-10,100)
    labels = ['Inner Radius', 'Outer Radius', 'Grain Size', 'Disk Mass', 'Surface Density', 'Beta']
    for i in range(6):
        plt.subplot(321+i)
        plt.plot(averages[i])
        plt.xlabel('Generation')
        plt.ylabel('Average ' + labels[i])
        plt.title(labels[i])
    """
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
                        wspace=None, hspace=None)
    """
    # DOES THIS OVERWRITE? IF NO, HANDLE DIFFERENTLY
    #plt.savefig('genetic/' + plotFile)
    
    plt.tight_layout()
    plt.show()
    
file_name = sys.argv[1]
process(file_name)
    from scipy import optimize
    import numpy as np
    import matplotlib.pyplot as plt
    import math, random
    from multiprocessing import Process, Pool

    #def fun(x):
    #    return [x[0]  + 0.5 * (x[0] - x[1])**3 - 1.0,0.5 * (x[1] - x[0])**3 + x[1]]
    #def jac(x):
    #    return np.array([[1 + 1.5 * (x[0] - x[1])**2,-1.5 * (x[0] - x[1])**2],[-1.5 * (x[1] - x[0])**2,1 + 1.5 * (x[1] - x[0])**2]])
    #sol = optimize.root(fun, [0, 0], jac=jac, method='hybr')
    #print(sol.x)

    K = 100      #Number of Biotic Components
    R = 100        #Essential Range (defines where Biotic Components can be present)
    P = 0          #Perturbation
    F = P
    OE = []        #Niche
    start = 0      #Time Start
    end = 200         #Time End
    step= 0.1    #Time Step
    w = []         #Affects Parameter (ranges between -1 and 1 for each K)
    u = []         #Ideal Temperature for species (between 0 and R -> the essential range)

    #OE = [random.uniform(3,10) for _ in range(K)] #Switches between same sized Niches to different sized ones
    OE = [5 for _ in range(K)]
    #populates affects values
    w = [random.uniform(-1,1) for _ in range(K)]

    while (int(len(w)) != int(len(set(w)))):
        print("Duplicate w's detected: Regenerating ...")
        w.clear()
        w = [random.uniform(-1,1) for _ in range(K)]

    #print(w)

    # #populates optimum growing temperatures [e.g 78.213423]
    u = [random.uniform(0, R) for _ in range(K)]

    while(int(len(u)) != int(len(set(u)))):
        print("Duplicate u's detected: Regenerating ...")
        u.clear()
        u = [random.uniform(0, R) for _ in range(K)]

    #print(u)

    N = 2           #Number of Environment Variables
    E = -20          #Temperature Start value
    Et = E         #Temperature without Biotic Force

    alpha = [[] for _ in range(K)] #abundance value for a species

    for _ in range(K):
        alpha[_].append((math.e) ** ((-1) * (((abs((E)-u[_])) ** 2) / (2*(OE[_]**2)))))

    rF = []         #Biotic Force Values
    rP = []         #Perturbation Values
    rE = []         #Temperature with Biotic Force Values
    rEt = []        #Temperature without Biotic Force Values

    #Abundance values over time
    rAx = [[] for x in range(K)]
    #Abundance values over time scaled up by R (Essential Range)
    rAxR = [[] for x in range(K)]
    #Tracks time (steps accumulation)
    time = []

    biotic_force = [[] for _ in range(K)]
    temperatures = []

    #plot abundance of species over temperature
    def plot_alphas():
        for x in np.arange (0, R, step):
            temperatures.append(x)

        for y in range(K):
            for x in np.arange (0, R, step):
                biotic_force[y].append((math.e) ** ((-1) * (((abs(x-u[y])) ** 2) / (2*(OE[y]**2)))) * w[y])

        plt.figure(figsize=(25,15))
        plt.title('Biotic Force over Time', fontsize=40)
        plt.xlabel('Temperature', fontsize=40)
        plt.ylabel('biotic force (a * w)', fontsize=40)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        for _ in range(K):
            plt.plot(temperatures,biotic_force[_])

        plt.plot(temperatures,np.sum((np.array(biotic_force, dtype=float)), axis=0), lw=4)
        plt.show()



    ########################################################################################################################
    def f1(x):
        biotic_force = []
        for y in range(K):
            biotic_force.append((math.e) ** ((-1) * (((abs(x-u[y])) ** 2) / (2*(OE[y]**2)))) * w[y])
        return(np.sum((np.array(biotic_force, dtype=float))))

    def stable_point_return(K_in):
        print("Running K", K_in)
        x = []
        y = []

        X1 = -50
        Y1 = R + 50

        for xi in np.arange(X1, Y1, 0.1):
            x.append(xi)
            y.append(f1(xi))

        true_zeros = []
        for _ in range(R):
            sol = optimize.root(f1, [_], jac=False, method='hybr')
            if(sol.x >=0 and sol.x <= R):
                true_zeros.append(sol.x)

        zeros_uniq = []
        for xi in true_zeros:
            if(int(xi) not in zeros_uniq):
                zeros_uniq.append(int(xi))

        print("Unique Points ...")
        zeros_uniq.sort()
        print(zeros_uniq)

        idx=0
        current_sign = "?"
        if(y[idx]>0):
            current_sign = "+"
        elif(y[idx]<0):
            current_sign = "-"

        stable_points = []

        for xi in np.arange(X1, Y1, 0.1):
            #print(x[idx])
            loopy_sign="?"
            #print(y[idx])
            if(y[idx]>0):
                loopy_sign = "+"
            elif(y[idx]<0):
                loopy_sign = "-"

            if(loopy_sign != current_sign):
                #print("Sign Change Detected!")
                #print("From : ", current_sign , " To : ", loopy_sign)
                if(current_sign == "+" and loopy_sign == "-"):
                    #print(">>>>> Stable Point : ", x[idx])
                    stable_points.append(int(x[idx]))
                current_sign = loopy_sign
                #print(y[idx])

            idx+=1


        print("Stable Points: ", stable_points)
        return(len(stable_points))

    biotic_elements_x = []
    stable_points_y   = []

    niche_vals = [5]

    stable_points_average = []

    def parallel_it(each_k):
        K = each_k
        ############################################################################
        w = [random.uniform(-1,1) for _ in range(K)]
        while (int(len(w)) != int(len(set(w)))):
            print("Duplicate w's detected: Regenerating ...")
            w.clear()
            w = [random.uniform(-1,1) for _ in range(K)]

        u = [random.uniform(0, R) for _ in range(K)]
        while(int(len(u)) != int(len(set(u)))):
            print("Duplicate u's detected: Regenerating ...")
            u.clear()
            u = [random.uniform(0, R) for _ in range(K)]

        alpha = [[] for _ in range(K)]
        for _ in range(K):
            alpha[_].append((math.e) ** ((-1) * (((abs((E)-u[_])) ** 2) / (2*(OE[_]**2)))))
        ############################################################################
        stable_points_average.append(stable_point_return(K))




    def plot_stable_biotic():
        plt.figure(figsize=(20,10))
        plt.title('Number of Species vs Stable Points', fontsize=40)
        plt.xlabel('Number of Species', fontsize=40)
        plt.ylabel('Number of Stable Points', fontsize=40)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.axvline(x=0)
        plt.axhline(y=0)
        plt.plot(biotic_elements_x,stable_points_y, 'r-',label = 'roots')
        #plt.savefig('stable_biotic.png')
        plt.show()


    def plot_function():
        print("Plotting Sum  ... ")
        plt.figure(figsize=(20,10))
        plt.title('xy', fontsize=40)
        plt.xlabel('x', fontsize=40)
        plt.ylabel('y', fontsize=40)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.axvline(x=0)
        plt.axhline(y=0)
        plt.plot(x,y, 'r-',label = 'roots')
        plt.show()


    def plot_stable_points():
        plt.figure(figsize=(20,10))
        plt.title('Stable Points', fontsize=40)
        plt.xlabel('temperature', fontsize=40)
        plt.ylabel('biotic force', fontsize=40)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        for stable in zeros_uniq:
            plt.axvline(x=stable)
        plt.axvline(x=0)
        plt.axhline(y=0)
        plt.plot(x,y, 'r-',label = 'temperature')
        plt.legend(loc=5, prop={'size': 30})
        plt.show()

    SAMPLE = 8
    if __name__ == '__main__':

        pool = Pool(processes=8)

        for each_niche in niche_vals:
            OE = [each_niche for _ in range(K)]
            for each_k in range(K):
                pool.map(parallel_it, [each_k for xp in range(SAMPLE)])

            biotic_elements_x.append(each_k)
            stable_points_y.append(sum(stable_points_average) / len(stable_points_average))
            stable_points_average.clear()
        #plot_alphas()
        #plot_function()
        #plot_stable_points()
        plot_stable_biotic()


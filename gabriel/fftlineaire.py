
def fft_gld_v2022(s,fe, fx):
    [dt, dy] = s.shape
    t = np.arange(0,dt)/fe
    Y = np.arange(0,dy)*fx
    debut = 1
    fin = len(Y)
    x = Y;
    Nmin = 1
    Nmax = dt
    print('Number of frames: ' + str(Nmax))
    H = np.zeros((Nmax-Nmin+1,fin-debut+1)); #initialisation du champ spatiotemporel total
    for i in range(Nmin,Nmax):
        H[i,:] = s[i,:]*fx
    compteur = 0; 
    Max_TF   = 10**(12);
    indicesx = np.zeros(len(x))

    ## TF temporelle 
    TFt = np.zeros((Nmax, len(x))); ## On initialise le tableau de la Transformee de fourier temporelle
    for j in range(1,len(x)):
        TFp = np.abs(fftshift(fft(H[:,j])-np.mean(H[i,:]))); ## Transformee de fourier spatiale 
        compteur = compteur + 1;
        indicesx[compteur]=j;
        TFt[:,compteur]=TFp;
    TFt[:,compteur] = np.abs(fftshift(fft(H[:,j])))
    ## Definition de l'axe des temps
    n = len(TFt[:,1])+1;
    f = np.arange(-fe/2,fe/2, fe/(n-1))

    midf=int(len(f)/2);
    # filtrage du mode n=1
    #TFt[midf-1:midf+1,:]=0;
    # 
    plt.figure()
    plt.plot(np.abs(f),TFt[:,1]**2);
    plt.xlabel('f')
    plt.ylabel('$TFt^2$')
    ## TF spatiale
    compteur = 0; 
    indicest = np.zeros(Nmax)
    TFx = np.zeros((Nmax, len(x))); # On initialise la transformee de fourier spatiale
    for i in range(1, Nmax):
        TFp=np.abs(fftshift(fft(H[i,:]-np.mean(H[i,:]))));
        compteur=compteur+1;
        indicest[compteur]=i;
        TFx[compteur,:]=TFp;

    #  Definition de l'axe des k
    ke = 1./fx; # echantillonage en k
    nx = len(TFx[1,:])+1;
    k  = np.arange(-ke/2,ke/2,ke/(nx-1))*2*np.pi; # axe des k 
    midk=int(len(k)/2);
    # filtrage du mode n=1
    #TFx[:,midk-1:midk+1]=0;
    # %% Moyenne temporelle de la transformee de fourier spatiale
    TF_x = np.mean(TFx, axis = 0);
    plt.figure()
    plt.loglog(np.abs(k),TF_x**2);
    plt.xlabel(r'k')
    plt.ylabel(r'$TFx^2$')
    
    TF=fftshift(fftshift(fft2(H-np.mean(H)),0),1); # fft dans l'axe des temps puis spatial

    nx=len(TF[1,:])+1;
    k=np.arange(-ke/2,ke/2,ke/(nx-1))*2*np.pi; # axe des k
    nt=len(TF[:,1])+1;
    f =np.arange(-fe/2,fe/2,fe/(nt-1)); # axe des temps
    midk=int(len(k)/2);
    midf=int(len(f)/2);
    # filtrage du mode n=1
    TF[:,midk-1:midk+1]=0;
    TF[midf-1:midf+1,:]=0; 
    TFT = np.abs(TF)
    
    # plt.figure()
    # plt.pcolormesh(k,f,TFT, vmin =0, vmax = 1000)
    # plt.xlabel(r'$k (\rm cm^{-1})$')
    # plt.ylabel(r'$f (\rm s^{-1})$')
    # plt.xlim([-5, 5])
    # plt.ylim([-25, 25])
    # cbar = plt.colorbar()
    # cbar.set_label(r'$S\left(\eta(y,t)\right)$')
    return k,f,TFT, TFt, TF_x
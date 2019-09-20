import numpy as np

def identity():
    return np.array([[1.,   0.],
                     [0.,   1.]])

def Mprop(phi):
    return np.array([[np.exp(1.j*phi), 0.],
                     [0., np.exp(-1.j*phi)]])

def Mbound(n1, n2):
    return 1./(2*n2) * np.array([[n2+n1, n2-n1],
                                 [n2-n1, n2+n1]])

def get_S_from_M(M):
    A = M[0,0]
    B = M[0,1]
    C = M[1,0]
    D = M[1,1]
    return 1./D * np.array([[A*D-B*C, B],
                            [-C,      1]])

def build_M(lst):
    M = identity()
    for i in lst:
        M = np.dot(i, M)
    return M

def get_list_of_matrices(wavelength, list_of_layers):

    k = 2.*np.pi/wavelength
    m = []   
    for i in np.arange(len(list_of_layers)):
        
        n = list_of_layers[i][0]
        d = list_of_layers[i][1]
        phi = n*k*d
        m.append(Mprop(phi))
        
        if i != len(list_of_layers)-1:
            n2 = list_of_layers[i+1][0]
            m.append(Mbound(n, n2))

    return m

def calc_t_and_r(M):
    S = get_S_from_M(M)
    return S[0,0], S[1,0]

def calc_t_and_r_from_layers(wavelength, list_of_layers):
    list_of_matrices = get_list_of_matrices(wavelength, list_of_layers)
    M = build_M(list_of_matrices)
    t,r = calc_t_and_r(M)
    return t,r

def get_field_pairs(wavelength, list_of_layers):

    num_layers = len(list_of_layers)    
    t,r = calc_t_and_r_from_layers(wavelength, list_of_layers)

    field_pairs = np.zeros((num_layers+1, 2)) +0.j

    E0 = 1.
    E0r = E0 * r
    field_pairs[0] = E0, E0r

    k = 2.*np.pi/wavelength

    for i in np.arange(num_layers):

        n1 = list_of_layers[i,0]
        phi = n1*k*list_of_layers[i,1]    

        if i != num_layers-1:

            n2 = list_of_layers[i+1,0]
            
            m = [Mprop(phi), Mbound(n1,n2)]
        else:
            m = [Mprop(phi)]
        M = build_M(m)

        field_pairs[i+1] = np.dot(M, field_pairs[i])

    return field_pairs

def get_field_in_between(wavelength, list_of_layers, field_pairs, N):

    k = 2.*np.pi/wavelength
    num_layers = len(list_of_layers)

    data_all_x = np.zeros((num_layers, N))
    data_all_y = np.zeros((num_layers, N)) +0.j

    distance = 0.
    
    for i in np.arange(num_layers):

        n, d = list_of_layers[i]

        x = np.linspace(0, d, N)
        y = field_pairs[i,0] * np.exp(1.j*n*k*x) + field_pairs[i,1] * np.exp(-1.j*n*k*x)

        data_all_x[i,:] = x + distance
        data_all_y[i,:] = y

        distance += d

    return data_all_x, data_all_y



def get_list_of_layers(n_SiO, d_SiO,
                       n_aSi, d_aSi,
                       num_pairs,
                       d_SiO_spacer,
                       d_aSi_spacer):

    list_of_layers = []

    for i in np.arange(num_pairs):
        list_of_layers.append([n_SiO, d_SiO])
        list_of_layers.append([n_aSi, d_aSi])

    list_of_layers.append([n_SiO, d_SiO_spacer])
    list_of_layers.append([n_aSi, d_aSi_spacer])
    list_of_layers.append([n_SiO, d_SiO_spacer])

    for i in np.arange(num_pairs):
        list_of_layers.append([n_aSi, d_aSi])
        list_of_layers.append([n_SiO, d_SiO])

    return np.array(list_of_layers)

def get_spectrum(list_of_layers, wavelengths):

    ts = np.copy(wavelengths)*0+0.j
##    rs = np.zeros((N))+0.j

    for i in np.arange(ts.size):
        wavelength = wavelengths[i]

        list_of_matrices = get_list_of_matrices(wavelength, list_of_layers)
        M = build_M(list_of_matrices)
        t, r = calc_t_and_r(M)

        ts[i] = t
##        rs[i] = r

    return np.abs(ts)**2
        



    
    


    

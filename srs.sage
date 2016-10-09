P.<x> = PolynomialRing(ZZ)
############################################

f = x^3-17*x^2-7*x-16 #sem napiste minimalny polynom (pokial mozno Pisotova cisla alebo aspon Salemova; inak to neskonci)
#f = x^5-45*x^4-39*x**3-33*x**2-3*x-27
#mp = 'plus' #'plus' alebo 'minus' beta rozvoje
#f=x^5-107*x^4-42*x^3 -42*x^2-42*x-42
mp = 'minus'

############################################
N.<beta> = NumberField(f, embedding = 100)
#f.complex_roots()
#map(abs,f.complex_roots())
#map(lambda x: x.N(), N.units())

def hasF(N,pm):
    '''
    overi, ze ma generator N vlastnost (plus/minus F), pm\in\{'plus','minus'\}
    '''
    b = N.gen()
#    print b.N()
    r = makevect(N,pm)
#    r = vector([0.22,0.2,-0.3,0.1])
    W = witnesses(r,N)
#    alpha = 0
    if pm == 'minus':
        alpha = b/(b+1)
#    alpha = 0.92
    print alpha.N()
    for w in W:
        print '-----------'
        out = SRS(w,r,N,alpha,1000)
        if out == 'perioda':
 #           continue
            return str(N.defining_polynomial()) + ' nema ' + pm + ' (F)'
        elif out == 'finite':
            continue
        else:
            return 'ani po 10000 iteracich to neskoncilo'
    return str(N.defining_polynomial()) + ' ma ' + pm + ' (F)'

############################################

def makevect(N,pm):
    '''
    vytvori vektor r prislusejici soustave
    '''
#    N.<beta> = NumberField(f, embedding = 100) #number field
    f = N.defining_polynomial()
    coef = f.coefficients(false)
    coef.reverse()
    d = f.degree()
    r = []
    for i in xrange(len(coef)):
        rs = vector(coef[i:],N)
        betas = vector([1/(beta^(k+1)) for k in xrange(len(rs))],N)
        r.append(rs*betas)
    if pm == 'plus':
        r = r[2:]
        r.reverse()
#        print map(lambda x: -x.N(),r)
        return -vector(r)
    r = [r[i]*(-1)^i for i in xrange(len(r))]
    r = r[2:]
    r.reverse()
    print map(lambda x: x.N(),r)
    print sum(r).N()
    return vector(r)
##########################################
def T(x,r,N,alpha):
    '''
    alpha-SRS transformace
    '''
#    beta = N.gen()
    new = list(x)
    new = new[1:]
    newentry = -(r*x+alpha).N().floor()
#    newentry = -(sum([r[i]*x[i] for i in xrange(len(r))])+alpha).N().floor()
    new.append(newentry)
    return vector(new)

def SRS(x,r,N,alpha,n):
    '''
    n-krat iterovane T, jestli najde netrivialni periodu, tak vyplivne 'perioda'; pokud se zacykli na nule, vyplivne 'finite'
    '''
    visited = [x]
    for i in xrange(n):
        print '-->' + str(x)
        if x == vector([0]*len(x)):
            return 'finite'
        x = T(x,r,N,alpha)
        
        if x in visited:
            print '--->' + str(visited[visited.index(x)]) + ' <---- uz bylo'
            return 'perioda'
        visited.append(x)
    return x
#######################################
def witnesses(r,N):
    '''
    konecna mnozina, kterou staci overit pro konecnost SRS
    '''
    dim = len(r)
    w = [[0]*i + [1] + [0]*(dim-i-1) for i in xrange(dim)]
    w = [vector(i) for i in w]
    W = w + [-i for i in w] #tady konecne mame +- standardni bazi
    i = 0
    while True:     #udelame tranzitivni uzaver na tau_{r,0}(z), -tau_{r,0}(-z)
        if i >= len(W):
            break
        new = [T(W[i],r,N,0)] + [-T(-W[i],r,N,0)]
        print str(W[i]) + ' -----> ' + str(new)
        if new[0] not in W:
            W.append(new[0])
        if new[1] not in W:
            W.append(new[1])
        i += 1
    print 'Velikost mnoziny svedku je ' + str(len(W))
    print 'Svedci jsou: ' + str(W)
    return W
makevect(N,'plus')
print hasF(N,mp)
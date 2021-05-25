#importando as bibliotecas necessarias 
import matplotlib.pyplot as plt        #biblioteca grafica
import numpy as np                     #biblioteca de operacoes matematicas
import scipy.optimize as optimization  #biblioteca de ajustes matematicos (least squares)

#inicio do programa

#definindo a funcao a ser ajustada:
#entrada:  x -> matriz que acomoda os vetores correspondentes a cada variavel
def func(x, b1, b2, b3, b4, b5, b6):
    return 4.*((b1*x[:,0] + b2*x[:,1] + b3*x[:,2] + b4*x[:,3])/(x[:,2] + b5*x[:,4] + 2.*b6*x[:,3]))

#funcao principal
def main():

    #abre e le os arquivos externos Pos e Nega
    arquivo1 = 'E_0_a_300_20_Co_20nm_B_Neg.txt'
    arquivo2 = 'E_0_a_300_20_Co_20nm_B_Pos.txt'
    a1 = np.loadtxt(arquivo1)
    a2 = np.loadtxt(arquivo2)

    t1 = np.zeros(len(a1))            #cria o vetor angulo 1
    t2 = np.radians(20.0)             #cria o vetor angulo 2
    intensidade1 = np.zeros(len(a1))  #cria o vetor intensidade arquivo 1
    intensidade2 = np.zeros(len(a2))  #cria o vetor intensidade arquivo 2
    dif_int = np.zeros(len(a1))       #cria o vetor diferenca de intensidade

    #armazena os arquivos externos nos vetores correspondentes
    for i in range(len(a1)):
        t1[i] = np.radians(a1[i,0])  #angulo em radianos
        intensidade1[i] = a1[i,1]    #intensidade arquivo 1
        intensidade2[i] = a2[i,1]    #intensidade arquivo2
        dif_int[i]=(intensidade2[i] - intensidade1[i])/intensidade2[i] #diferenca de intensidade

    f1 = np.zeros(len(t1)) #cria vetor para a funcao f1
    f2 = np.zeros(len(t1)) #cria vetor para a funcao f2
    f3 = np.zeros(len(t1)) #cria vetor para a funcao f3
    f4 = np.zeros(len(t1)) #cria vetor para a funcao f4
    f5 = np.zeros(len(t1)) #cria vetor para a funcao f5

    #loop para atribuir valor a cada funcao fi
    for i in range(len(t1)):
        f1[i] = ((np.sin(t1[i]))**2)*np.sin(t2)*np.cos(t2) - ((np.sin(t2))**2)*np.sin(t1[i])*np.cos(t1[i]) 
        f2[i] = ((np.cos(t2))**2)*np.sin(t1[i])*np.cos(t1[i]) - ((np.cos(t1[i]))**2)*np.sin(t2)*np.cos(t2)
        f3[i] = ((np.sin(t1[i]))**2)*((np.sin(t2))**2)
        f4[i] = np.sin(t1[i])*np.cos(t1[i])*np.sin(t2)*np.cos(t2)
        f5[i] = ((np.cos(t1[i]))**2)*((np.cos(t2))**2)

    #armazena as funcoes fi em uma matriz (nx6)
    m = np.transpose(np.matrix([f1,f2,f3,f4,f5]))
    
    #x0 = [1.0,1.0,1.0,1.0,1.0,1.0]
    
    #executa o metodo de ajuste por minimos quadrados
    #entra: func    --> nome da funcao a ser ajustada
    #       m       --> matriz com as variaveis independentes (x1,x2,x3,...)
    #       dif_int --> vetor da variavel dependente (y)
    #sai:   c1 --> matriz com os melhores coeficientes
    #       c2 --> matriz de covariancia
    c1,c2=optimization.curve_fit(func, m, dif_int)
    
    #armazena os valores das constantes ajustadas
    b1=c1[0]
    b2=c1[1]
    b3=c1[2]
    b4=c1[3]
    b5=c1[4]
    b6=c1[5]

    #cria uma funcao para visualizacao grafica do ajuste feito
    y = np.zeros(len(t1))
    for i in range(len(t1)):
        y[i] = 4.*((b1*m[i,0] + b2*m[i,1] + b3*m[i,2] + b4*m[i,3])/(m[i,2] + b5*m[i,4] + 2.*b6*m[i,3]))

    #plota os dados nao tratados e o ajuste feito.
    plt.plot(t1,dif_int,'p', label='Dados')
    plt.plot(t1,y,color='red', label='Ajuste')
    plt.legend()
    plt.show()

    return

main()

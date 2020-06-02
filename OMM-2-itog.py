import numpy as np
import matplotlib.pyplot as plt

Nx, Ny, Nt = 50, 50, 50 # количество шагов по x, y и t
Xmax = np.pi/2 # максимальное значение x
Ymax = np.pi # максимальное значение y
Tmax = 1  # максимальное значение t
          # - будем менять, чтобы посмотреть как меняется при этом график

X0, Y0, T0 = 0, 0, 0 # начальные значения x,y и t соответственно
hx = (Xmax-X0)/Nx # шаг по x
hy = (Ymax-Y0)/Ny # шаг по y
ht = (Tmax-T0)/Nt # шаг по t

# заполняем сетку
x = np.linspace(X0, Xmax, Nx)
y = np.linspace(Y0, Ymax, Ny)
t = np.linspace (T0, Tmax, Nt)

 # заполняем решение нулевыми значениями
w = np.zeros((Nx, Ny, 2 * Nt + 1))  # массив по t имеет размер 2Nt+1 так как при каждом переходе с слоя на слой мы добавляем промежутчный слой k+1/2

# создаем функции F_1 и F_2 - в файле описания это функции F(j+1/2) и F(j+1)
def F_1(i, k, j):
    return (9*ht/(2*hy**2))*(w[i, k+1, j-1] + w[i, k-1, j-1])-((9*ht/(hy**2))-1)*w[i, k, j-1]

def F_2(i, k, j):
    return (9*ht/(2*hx**2))*(w[i+1, k, j-1] + w[i-1, k, j-1])-((9*ht/(hx**2))-1)*w[i, k, j-1]

# прогонка
def progonka_x (k, j):
    # Вводим массивы для прогоночных коэффициентов α и β
    α = np.zeros(Nx-1)
    β = np.zeros(Nx-1)
    α[0] = 0 # из граничных условий
    β[0] = 0
    A = (9*ht)/(2*hx**2)
    B = (9*ht)/(2*hx**2)
    C = 1 + 9*ht/(hx**2)
    # Заполняем массив коэффициентов α и β по формулам для  β и α
    for l in range(1, Nx-2):
        Fl = -F_1(l, k, j)
        α[l+1] = B/(C - A * α[l])
        β[l+1] = (Fl + A * β[l])/(C - A * α[l])
    # из граничного условия yN = 0
    w[Nx-1, k, j] = β[-1]/(1 - α[-1])
    for l in range(Nx-2, 0, -1):
        w[l-1, k, j] = α[l]*w[l, k, j] + β[l]

def progonka_y (i, j):
    α = np.zeros(Ny - 1)
    β = np.zeros(Ny - 1)
    α[0] = 0
    β[0] = 0
    A = (9*ht)/(2*hy**2)
    B = (9*ht)/(2*hy**2)
    C = 1 + 9*ht/(hy**2)
    for l in range(1, Ny-2):
        Fl = -F_2(i, l, j)
        α[l+1] = B/(C - A * α[l])
        β[l+1] = (Fl + A * β[l])/(C - A * α[l])
    # из граничного условия yN = 0
    w[i, Ny-1, j] = β[-1]/(1 - α[-1])
    for l in range(Ny-2, 0, -1):
        w[i, l-1, j] = α[l]*w[i, l, j] + β[l]

w[:, :, 0] = np.sin(2*x)*np.sin(3*y) #начальное условие


for j in range(1, 2*Nt, 2):
    for k in range(1, Ny - 1):
        progonka_x(k, j) # переходим на слой 1/2
    for i in range(1, Nx - 1):
        progonka_y(i, j+1)

# массив аналитического решения
u = np.zeros((Nx, Ny))
for k in range(1, Nx):
    for i in range(1, Ny):
        u[i,k] = np.exp(-117*t[1])*np.sin(2*x[i])*np.sin(3*y[k])

# строим графики  Численного решения и Аналитического решений

y, x = np.meshgrid(y, x)

fig = plt.figure()
chislennoe = fig.add_subplot(121, projection='3d')
analit = fig.add_subplot(122, projection='3d')
surf1 = chislennoe.plot_surface(x, y, w[:, :, 1])
surf2 = analit.plot_surface(x, y, u[:, :])

# С названиями осей пока не заморачивалась, слева численное решение, а справа - аналитическое. В идеале должны быть одинаковые

plt.show()

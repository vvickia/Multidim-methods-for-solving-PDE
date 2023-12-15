# Многомерные методы решения уравнений в частных производных

Программа, реализующая численное решение двумерного уравнения теплопроводности с помощью _полностью неявной схемы_.

Для расщепления двумерной схемы по пространственным направлениям выбрана продольно-поперечная схема (метод переменных направлений).

## Постановка задачи

Промоделировать с помощью разработанной программы установление температуры в однородной прямоугольной медной пластине размерами $L_x \times L_y = 0.5 \times 0.5$ м. Торцы пластины вдоль одной из осей поддерживаются при постоянной температуре $25^\circ C$. На торцах вдоль другой оси задан нулевой поток тепла. Начальная температуры пластины равна $0^\circ C$.

Построить рисунки с двумерными картами распределения температуры для нескольких выбранных моментов времени. Проанализировать характер решения.
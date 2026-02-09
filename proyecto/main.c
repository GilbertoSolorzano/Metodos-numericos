/*
Nombre del programa: Proyecto final Metodos Numericos
Autor: Solorzano Galvez Gilberto Jesus
Fecha: Mayo 2024
Hardware: Core i5-6300u 2.50 GHz
Objetivo: este programa cuenta con los diferentes metodos de los metodos numericos que se vieron a lo largo del semestre,
ejecuta las diferentes funciones pidiendo en algunos casos los valores de x
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

double fun(double x){
	double y=2*x*x*x - 9*x*x + 12*x - 3;
	return y;
}
double funDev(double x){
	double y=sin(x);
	return y;	
}
double fpun(double x){
	double y = x*x - 0.5;
	return y;
}
double fun2(double x){
	double y=6*x*x - 18*x + 12;
	return y;
}
double fun3(double x){
	double y = 12*x -18;
	return y;
}
double f(double x, double y){
	double z = sin(x) + y;
	return z;
}
void tanteos(){
    double delta, e, x, y;
    printf("Ingrese el valor del salto: ");
    scanf("%lf", &delta);
    printf("Ingrese el valor del error: ");
    scanf("%lf", &e);

    x = delta;
    for(int i = 0; i < 10000; i++){
        y = fun(x);
        if(y == 0){
            printf("%lf es un aproximado a la raiz", x);
            return;
        }else{
            if(y < 0){
                if(delta <= e){
                    printf("%lf es un aproximado a la raiz", x);
                    return;
                }else{
                    x = x - delta;
                    delta = delta / 10;
                }
            }
            x = x + delta;
        }
    }
}
void punto_fijo() {
    
    
    double nmi, e, x, y, rel;
    printf("Ingresa el valor de X:");
    scanf("%lf", &x);
    printf("Ingresa el valor del error:");
    scanf("%lf", &e);
    printf("Ingresa el numero de iteraciones:");
    scanf("%lf", &nmi);
	
	clock_t start, end;
    start = clock();
	
    for (int i = 0; i <= nmi; i++) {
    	
        printf("entro %lf \n", x);
        y = fpun(x) + x;
        rel = (x - y) / y;
        if (rel < 0)
            rel *= -1;
        printf("y= %lf \nrel= %lf\n", y, rel);
        if (rel <= e) {
            end = clock();
            double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
            printf("Tiempo de ejecucion: %f segundos\n", time_taken);
            printf("Si converge en %lf iteraciones, y= %lf\n", nmi, y);
            return;
        } else {
            x = y;
        }
    }
    
    end = clock();
    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Tiempo de ejecucion: %f segundos\n", time_taken);
    printf("No converge en %lf iteraciones\n", nmi);
    return;
}
void nr_1er() {
    int nmi;
    double y, yp, rel, x, e;
    
    printf("Ingresa el valor de X: ");
    scanf("%lf", &x);
    printf("Ingresa el valor del error: ");
    scanf("%lf", &e);
    printf("Ingresa el numero de iteraciones: ");
    scanf("%d", &nmi);

    clock_t start, end;
    start = clock();

    for (int i = 0; i <= nmi; i++) {
        y = fun(x);
        yp = fun2(x);
        y = x - (y / yp);
        rel = (y - x) / y;
        if (rel < 0)
            rel *= -1;
        if (rel <= e) {
            end = clock();
            double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
            printf("Tiempo de ejecucion: %f segundos\n", time_taken);
            printf("Una aproximacion a la raiz %lf\n", y);
            return;
        } else {
            x = y;
        }
    }

    end = clock();
    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Tiempo de ejecucion: %f segundos\n", time_taken);
    printf("No converge en %d calculos\n", nmi);
    return;
}
void newton_2do() {
    int delta, nmi;
    double e, x, y, yp, rel;
    
    printf("Ingresa el valor del X: ");
    scanf("%lf", &x);
    printf("Ingresa el valor del error: ");
    scanf("%lf", &e);
    printf("Ingresa el numero de iteraciones: ");
    scanf("%d", &nmi);

    clock_t start, end;
    start = clock();

    for (int i = 0; i <= nmi; i++) {
        y = x - ((fun(x) * fun2(x)) / ((fun2(x) * fun2(x)) - (fun(x) * fun3(x))));
        rel = (y - x) / y;
        printf("x= %lf\n y= %lf\n rel=%lf\n", x, y, rel);
        if (rel < 0)
            rel *= -1;
        if (rel <= e) {
            end = clock();
            double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
            printf("Tiempo de ejecucion: %f segundos\n", time_taken);
            printf("Una aproximacion a la raiz %lf\n", y);
            return;
        } else {
            x = y;
        }
    }

    end = clock();
    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Tiempo de ejecucion: %f segundos\n", time_taken);
    printf("No converge en %d iteraciones\n", nmi);
    return;
}
void biseccion() {
    double x0, x1, e, y0, y1, xm, ym, delta;
    int i = 0;

    printf("Ingresa el valor de X0: ");
    scanf("%lf", &x0);
    printf("Ingresa el valor de X1: ");
    scanf("%lf", &x1);
    printf("Ingresa el valor del error: ");
    scanf("%lf", &e);

    y0 = fun(x0);
    y1 = fun(x1);

    if (y0 * y1 >= 0) {
        printf("Para que funcione, el producto de la multiplicacion de x0 y x1 debe ser menor a 0, cambiar los valores\n");
        biseccion();
    } else {
        clock_t start, end;
        start = clock();

        do {
            xm = (x0 + x1) / 2;
            ym = fun(xm);

            printf("\nIteracion %d:\n", i);
            printf("x0 = %lf\n", x0);
            printf("x1 = %lf\n", x1);
            printf("xm = (x0 + x1)/2 = %lf\n", xm);
            printf("ym = fun(xm) = %lf\n", ym);

            if (y0 * ym < 0) {
                x1 = xm;
                printf("Cambiando x1 = xm\n");
            } else {
                x0 = xm;
                y0 = ym;
                printf("Cambiando x0 = xm\n");
            }

            delta = (x1 - x0) / x0;
            if (delta < 0)
                delta *= -1;
            
            printf("delta = (x1 - x0)/x0 = %lf\n", delta);

            i++;
        } while (!(delta <= e));

        end = clock();
        double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
        printf("Tiempo de ejecucion: %f segundos\n", time_taken);
        printf("\nUn aproximado a la raiz es: %lf\n", xm);
    }
}
void falsa_pos() {
    double x0, x1, e, y0, y1, xm, ym, delta;
    
    printf("Ingresa el valor de X0: ");
    scanf("%lf", &x0);
    printf("Ingresa el valor de X1: ");
    scanf("%lf", &x1);
    printf("Ingresa el valor del error: ");
    scanf("%lf", &e);
    
    y0 = fun(x0);
    y1 = fun(x1);
    
    if (y0 * y1 >= 0) {
        printf("Para que funcione, el producto de la multiplicacion de x0 y x1 debe ser menor a 0, cambiar los valores\n");
        falsa_pos();
    } else {
        clock_t start, end;
        start = clock();

        do {
            y0 = fun(x0);
            y1 = fun(x1);
            xm = x1 - ((y1 * (x1 - x0)) / (y1 - y0));
            ym = fun(xm);
            delta = (xm - x0) / xm;
            if (delta < 0)
                delta *= -1;
            if (y0 * ym < 0) {
                x1 = xm;
                y1 = ym;
            } else {
                x0 = xm;
                y0 = ym;
            }
        } while (delta >= e);

        end = clock();
        double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
        printf("Tiempo de ejecucion: %f segundos\n", time_taken);
        printf("Un aproximado a la raiz es: %lf\n", xm);
    }
    return;
}
void secante() {
    double x0, x1, e, nmi, y0, y1, xm, delta;
    
    printf("Ingresa el valor de X0: ");
    scanf("%lf", &x0);
    printf("Ingresa el valor de X1: ");
    scanf("%lf", &x1);
    printf("Ingresa el valor del error: ");
    scanf("%lf", &e);
    printf("Ingresa el numero de iteraciones: ");
    scanf("%lf", &nmi);

    clock_t start, end;
    start = clock();

    for (int i = 0; i <= nmi; i++) {
        y0 = fun(x0);
        y1 = fun(x1);
        xm = x1 - ((y1 * (x0 - x1)) / (y0 - y1));
        delta = (xm - x1) / xm;
        if (delta < 0)
            delta *= -1;
        if (delta <= e) {
            end = clock();
            double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
            printf("Tiempo de ejecucion: %f segundos\n", time_taken);
            printf("Un aproximado a la raiz es: %lf\n", xm);
            return;
        } else {
            x0 = x1;
            x1 = xm;
        }
    }

    end = clock();
    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Tiempo de ejecucion: %f segundos\n", time_taken);
    printf("No converge en %lf calculos\n", nmi);
    return;
}

void jacobi() {
    double x1 = 0, x2 = 0, x3 = 0, nmi, delta1, delta2, delta3, e, i = 0, y1, y2, y3;
    
    printf("Ingresa el valor del error: ");
    scanf("%lf", &e);
    printf("Ingresa el numero de iteraciones: ");
    scanf("%lf", &nmi);

    clock_t start, end;
    start = clock();

    do {
        y1 = (20 + x2 - x3) / 20;
        y2 = (11 - 2 * x1 + x3) / 10;
        y3 = (-18 - x1 - x2) / (-20);
        
        delta1 = fabs((y1 - x1) / y1);
        delta2 = fabs((y2 - x2) / y2);
        delta3 = fabs((y3 - x3) / y3);

        x1 = y1;
        x2 = y2;
        x3 = y3;
        
        i++;

        if (delta1 <= e && delta2 <= e && delta3 <= e) {
            end = clock();
            double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
            printf("Tiempo de ejecucion: %f segundos\n", time_taken);
            printf("Aproximacion a las raices:\n y1 = %lf \n y2 = %lf \n y3 = %lf \n", y1, y2, y3);
            return;
        }
        
    } while (i < nmi);

    end = clock();
    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Tiempo de ejecucion: %f segundos\n", time_taken);
    printf("No converge\n");
}
void seidell() {
    double x1 = 0, x2 = 0, x3 = 0, nmi, delta, e, i = 0, x0 = 0;
    
    printf("Ingresa el valor del error: ");
    scanf("%lf", &e);
    printf("Ingresa el numero de iteraciones: ");
    scanf("%lf", &nmi);

    clock_t start, end;
    start = clock();

    do {
        x1 = (20 + x2 - x3) / 20;
        x2 = (11 - 2 * x1 + x3) / 10;
        x3 = (-18 - x1 - x2) / (-20);
        delta = (x1 - x0) / x1;
        if (delta < 0)
            delta *= -1;
        x0 = x1;
        i++;
        if (delta <= e) {
            end = clock();
            double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
            printf("Tiempo de ejecucion: %f segundos\n", time_taken);
            printf("Aproximacion a las raices:\n x1 = %lf \n x2 = %lf \n x3 = %lf \n", x1, x2, x3);
            return;
        }
        
    } while (i < nmi);

    end = clock();
    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Tiempo de ejecucion: %f segundos\n", time_taken);
    printf("No converge\n");
}
//diferenciacion 
void pDerivdabaja() {
    double h, x, x1, x2, res; 

    printf("Ingresa el valor de X: ");
    scanf("%lf", &x); 
    printf("Ingresa el valor de h (salto): ");
    scanf("%lf", &h); 

    if (h == 0) {
        printf("El valor de h no puede ser cero.\n");
        return;
    }

    x1 = x + h;
    x2 = x - h;

    clock_t start, end;
    start = clock();

    res = (funDev(x1) - funDev(x2)) / (2 * h); // Cambio de funcion por funDev
    
    end = clock();
    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Tiempo de ejecucion: %f segundos\n", time_taken);

    printf("El resultado de la aproximacion de la derivada es: %lf\n", res);
} 
void pDerivda() {
    double h, x, res;

    printf("Ingresa el valor de X: ");
    scanf("%lf", &x);
    printf("Ingresa el valor de h (salto): ");
    scanf("%lf", &h);

    if (h == 0) {
        printf("El valor de h no puede ser cero.\n");
        return;
    }

    clock_t start, end;
    start = clock();

    double x1 = x - (2 * h);
    double x2 = x - h;
    double x3 = x;
    double x4 = x + h;
    double x5 = x + (2 * h);

    printf("\tValores de X\n");
    printf("Xi-2= %lf\n    Xi-1= %lf\n    Xi= %lf\n  Xi+1= %lf\n   Xi+2= %lf\n\n", x1, x2, x3, x4, x5);

    double a = funDev(x1);
    double b = funDev(x2);
    double c = funDev(x3);
    double d = funDev(x4);
    double e = funDev(x5);

    printf("\tValores f(x)\n");
    printf("f(Xi-2) = %lf\n f(Xi-1)=%lf\n f(Xi)=%lf\n f(Xi+1)=%lf\n f(Xi+2)=%lf\n\n", a, b, c, d, e);

    res = ((-1 * funDev(x5)) + 8 * funDev(x4) - 8 * funDev(x2) + funDev(x1)) / (12 * h);

    end = clock();
    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Tiempo de ejecucion: %f segundos\n", time_taken);

    printf("El resultado de la aproximacion de la derivada es: %lf\n", res);
}
void sDerivda() {
    double h, x, res;

    printf("Ingresa el valor de X: ");
    scanf("%lf", &x);
    printf("Ingresa el valor de h (salto): ");
    scanf("%lf", &h);

    if (h == 0) {
        printf("El valor de h no puede ser cero.\n");
        return;
    }

    clock_t start, end;
    start = clock();

    double x1 = x - (2 * h);
    double x2 = x - h;
    double x3 = x;
    double x4 = x + h;
    double x5 = x + (2 * h);
    res = (-funDev(x5) + 16 * funDev(x4) - 30 * funDev(x3) + 16 * funDev(x2) - funDev(x1)) / (12 * h * h);

    end = clock();
    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Tiempo de ejecucion: %f segundos\n", time_taken);

    printf("El resultado de la aproximación de la segunda derivada es: %lf\n", res); 
}
void tDerivda() {
    double h, x, res;

    printf("Ingresa el valor de X: ");
    scanf("%lf", &x);
    printf("Ingresa el valor de h (salto): ");
    scanf("%lf", &h);

    if (h == 0) {
        printf("El valor de h no puede ser cero.\n");
        return;
    }

    clock_t start, end;
    start = clock();

    double x1 = x - (3 * h);
    double x2 = x - (2 * h);
    double x3 = x - h;
    double x4 = x;
    double x5 = x + h;
    double x6 = x + (2 * h);
    double x7 = x + (3 * h);
    res = (-(funDev(x7)) + 8 * (funDev(x6)) - 13 * (funDev(x5)) + 13 * (funDev(x3)) - 8 * (funDev(x2)) + funDev(x1)) / (8 * (h * h * h));

    end = clock();
    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Tiempo de ejecucion: %f segundos\n", time_taken);

    printf("El resultado de la aproximacion de la tercera derivada es: %lf\n", res); 
}
void cDerivada() {
    double h, x, res;

    printf("Ingresa el valor de X: ");
    scanf("%lf", &x);
    printf("Ingresa el valor de h (salto): ");
    scanf("%lf", &h);

    if (h == 0) {
        printf("El valor de h no puede ser cero.\n");
        return;
    }

    clock_t start, end;
    start = clock();

    double x1 = x - (3 * h);
    double x2 = x - (2 * h);
    double x3 = x - h;
    double x4 = x;
    double x5 = x + h;
    double x6 = x + (2 * h);
    double x7 = x + (3 * h);
    res = (-(funDev(x7)) + 12 * (funDev(x6)) - 39 * (funDev(x5)) + 56 * funDev(x4) - 39 * (funDev(x3)) + 12 * (funDev(x2)) - funDev(x1)) / (6 * (h * h * h * h));

    end = clock();
    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Tiempo de ejecucion: %f segundos\n", time_taken);

    printf("El resultado de la aproximacion de la cuarta derivada es: %lf\n", res); 
}

//integrales
void trapeSimple(){
    double b, a, x0,x1, y;
    printf("Ingrese el intervalo cerrado para la integral\n a: ");
    scanf("%lf", &a);
    printf("b: ");
    scanf("%lf", &b);
    x0=a;
    x1=b;

    clock_t start, end;
    start = clock();

    y = ((b-a)/2)*(fun(x0)+fun(x1));

    end = clock();
    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Tiempo de ejecucion: %f segundos\n", time_taken);

    printf("el resultado de la integral es: %lf \n", y);
}
void trapemulti(){
    double h, x0, xn, xi, n, sum = 0, a, b, y;
    printf("Ingrese el intervalo cerrado para la integral\n a: ");
    scanf("%lf", &a);
    printf("b: ");
    scanf("%lf", &b);
    printf("Ingrese el numero de segmentos (2,3,4,...):");
    scanf("%lf", &n);
    
    x0=a;
    xi=x0;
    xn=b;
    h=(b-a)/n;
    
    clock_t start, end;
    start = clock();

    for(int i=1; i<=n-1; i++){
        sum += fun(xi + h*i);
    }
    y=(h/2)*(fun(x0)+fun(xn)+2*sum);

    end = clock();
    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Tiempo de ejecucion: %f segundos\n", time_taken);

    printf("el resultado de la integral es: %lf \n", y);
}
void simpson13Simple(){
    double b, a, x0,x1,x2, y, h;
    printf("Ingrese el intervalo cerrado para la integral\n a: ");
    scanf("%lf", &a);
    printf("b: ");
    scanf("%lf", &b);
    h=(b-a)/2;
    x0=a;
    x1=x0+h;
    x2=x1+h;

    clock_t start, end;
    start = clock();

    y=((b-a)/6)*(fun(x0) + 4*fun(x1) + fun(x2));    

    end = clock();
    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Tiempo de ejecucion: %f segundos\n", time_taken);

    printf("el resultado de la integral es: %lf\n", y);
}
void simpson13multi() {
    double h, x0, xn, xi, n, sum1 = 0, sum2 = 0, a, b, y;
    int i;

    printf("Ingrese el intervalo cerrado para la integral\n a: ");
    scanf("%lf", &a);
    printf("b: ");
    scanf("%lf", &b);
    printf("Ingrese el numero de segmentos (debe ser par): ");
    scanf("%lf", &n);

    // Verificar que n es par
    if ((int)n % 2 != 0) {
        printf("El numero de segmentos debe ser par.\n");
        return;
    }

    x0 = a;
    xn = b;
    h = (b - a) / n;

    clock_t start, end;
    start = clock();

    // Sumar las funciones evaluadas en los puntos impares
    for (i = 1; i <= n - 1; i += 2) {
        sum1 += fun(x0 + i * h);
    }

    // Sumar las funciones evaluadas en los puntos pares
    for (i = 2; i <= n - 2; i += 2) {
        sum2 += fun(x0 + i * h);
    }

    y = (h / 3) * (fun(x0) + fun(xn) + 4 * sum1 + 2 * sum2);

    end = clock();
    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Tiempo de ejecucion: %f segundos\n", time_taken);

    printf("El resultado de la integral es: %lf\n", y);
}
void simpson38Simple(){
    double b, a, x0,x1,x2,x3, y, h;
    printf("Ingrese el intervalo cerrado para la integral\n a: ");
    scanf("%lf", &a);
    printf("b: ");
    scanf("%lf", &b);
    h=(b-a)/3;
    x0=a;
    x1=x0+h;
    x2=x1+h;
    x3=x2+h;

    clock_t start, end;
    start = clock();

    y=((b-a)/8)*(fun(x0)+fun(x3) + 3*fun(x1) + 3*fun(x2));    

    end = clock();
    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Tiempo de ejecucion: %f segundos\n", time_taken);

    printf("El resultado de la integral es: %lf\n", y);    
}
void simpson38multi() {
    double h, x0, xn, xi, n, sum1 = 0, sum2 = 0, sum3 = 0, a, b, y;
    int i;

    printf("Ingrese el intervalo cerrado para la integral\n a: ");
    scanf("%lf", &a);
    printf("b: ");
    scanf("%lf", &b);
    printf("Ingrese el numero de segmentos (multiplo de 3): ");
    scanf("%lf", &n);

    // Verificar que n es múltiplo de 3
    if ((int)n % 3 != 0) {
        printf("El número de segmentos debe ser multiplo de 3.\n");
        return;
    }

    x0 = a;
    xn = b;
    h = (b - a) / n;

    clock_t start, end;
    start = clock();

    // Sumar las funciones evaluadas en los puntos intermedios
    for (i = 1; i <= n - 1; i++) {
        xi = x0 + i * h;
        if (i % 3 == 0) {
            sum3 += fun(xi);
        } else {
            sum1 += fun(xi);
        }
    }

    y = (3 * h / 8) * (fun(x0) + fun(xn) + 2 * sum3 + 3 * sum1);

    end = clock();
    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Tiempo de ejecucion: %f segundos\n", time_taken);

    printf("El resultado de la integral es: %lf\n", y);
}
void boole245Simple() {
    double b, a, x0, x1, x2, x3, x4, y, h;
    printf("Ingrese el intervalo cerrado para la integral\n a: ");
    scanf("%lf", &a);
    printf("b: ");
    scanf("%lf", &b);
    h = (b - a) / 4;
    x0 = a;
    x1 = x0 + h;
    x2 = x1 + h;
    x3 = x2 + h;
    x4 = x3 + h;

    clock_t start, end;
    start = clock();

    y = ((b - a) / 90) * (7 * fun(x0) + 7 * fun(x4) + 32 * fun(x3) + 32 * fun(x1) + 12 * fun(x2));

    end = clock();
    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Tiempo de ejecucion: %f segundos\n", time_taken);

    printf("El resultado de la integral es: %lf\n", y);
}
void boole245multi() {
    double h, x0, xn, xi, n, sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, a, b, y;
    int i;

    printf("Ingrese el intervalo cerrado para la integral\n a: ");
    scanf("%lf", &a);
    printf("b: ");
    scanf("%lf", &b);
    printf("Ingrese el numero de segmentos (multiplo de 4): ");
    scanf("%lf", &n);

    // Verificar que n es múltiplo de 4
    if ((int)n % 4 != 0) {
        printf("El número de segmentos debe ser multiplo de 4.\n");
        return;
    }

    x0 = a;
    xn = b;
    h = (b - a) / n;

    clock_t start, end;
    start = clock();

    // Sumar las funciones evaluadas en los puntos intermedios
    for (i = 1; i <= n - 1; i++) {
        xi = x0 + i * h;
        if (i % 4 == 1) {
            sum1 += fun(xi);
        } else if (i % 4 == 2) {
            sum2 += fun(xi);
        } else if (i % 4 == 3) {
            sum3 += fun(xi);
        } else {
            sum4 += fun(xi);
        }
    }

    y = (2 * h / 45) * (7 * fun(x0) + 7 * fun(xn) + 32 * sum1 + 12 * sum2 + 32 * sum3 + 14 * sum4);

    end = clock();
    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Tiempo de ejecucion: %f segundos\n", time_taken);

    printf("El resultado de la integral es: %lf\n", y);
}
void boole5288Simple() {
    double b, a, x0, x1, x2, x3, x4, x5, y, h;
    printf("Ingrese el intervalo cerrado para la integral\n a: ");
    scanf("%lf", &a);
    printf("b: ");
    scanf("%lf", &b);
    h = (b - a) / 5;
    x0 = a;
    x1 = x0 + h;
    x2 = x1 + h;
    x3 = x2 + h;
    x4 = x3 + h;
    x5 = x4 + h;

    clock_t start, end;
    start = clock();

    y = ((b - a) / 288) * (19 * (fun(x0) + fun(x5)) + 75 * (fun(x1) + fun(x4)) + 50 * (fun(x2) + fun(x3)));

    end = clock();
    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Tiempo de ejecucion: %f segundos\n", time_taken);

    printf("El resultado de la integral es: %lf\n", y);
}
void boole5288multi() {
    double h, x0, xn, xi, n, sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, sum5 = 0, a, b, y;
    int i;

    printf("Ingrese el intervalo cerrado para la integral\n a: ");
    scanf("%lf", &a);
    printf("b: ");
    scanf("%lf", &b);
    printf("Ingrese el numero de segmentos (multiplo de 5): ");
    scanf("%lf", &n);

    // Verificar que n es múltiplo de 5
    if ((int)n % 5 != 0) {
        printf("El numero de segmentos debe ser multiplo de 5.\n");
        return;
    }

    x0 = a;
    xn = b;
    h = (b - a) / n;

    clock_t start, end;
    start = clock();

    // Sumar las funciones evaluadas en los puntos intermedios
    for (i = 1; i <= n - 1; i++) {
        xi = x0 + i * h;
        if (i % 5 == 1) {
            sum1 += fun(xi);
        } else if (i % 5 == 2) {
            sum2 += fun(xi);
        } else if (i % 5 == 3) {
            sum3 += fun(xi);
        } else if (i % 5 == 4) {
            sum4 += fun(xi);
        } else {
            sum5 += fun(xi);
        }
    }

    y = (5 * h / 288) * (19 * (fun(x0) + fun(xn)) + 75 * (sum1 + sum5) + 50 * (sum2 + sum4) + 38 * sum3);

    end = clock();
    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Tiempo de ejecucion: %f segundos\n", time_taken);

    printf("El resultado de la integral es: %lf\n", y);
}
//ecuaciones diferenciales
void euler() {
    double x0, xf, y, h, n;
    printf("Ingresa el valor de X0: ");
    scanf("%lf", &x0);
    
    printf("Ingresa el valor de XF: ");
    scanf("%lf", &xf);

    printf("Ingresa el valor de Y: ");
    scanf("%lf", &y);
    
    printf("Ingresa el numero de iteraciones: ");
    scanf("%lf", &n);

    h = (xf - x0) / n;

    clock_t start, end;
    start = clock();

    // Implementación del método de Euler
    for (int i = 0; i < n; ++i) {
        y = y + h * f(x0, y);
        x0 = x0 + h;
    }

    end = clock();
    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Tiempo de ejecucion: %f segundos\n", time_taken);

    printf("El valor de Y en XF es: %lf\n", y);
}
void eulerGauss() {
    double x, y, yp, yc, rel, h, n, a, b, e;
    clock_t start, end;
    double cpu_time_used;

    printf("Ingresa el valor de X: ");
    scanf("%lf", &x);

    printf("Ingresa el valor de a: ");
    scanf("%lf", &a);

    printf("Ingresa el valor de b: ");
    scanf("%lf", &b);

    printf("Ingresa el valor de Y: ");
    scanf("%lf", &y);

    printf("Ingresa el numero de iteraciones: ");
    scanf("%lf", &n);

    printf("Ingresa el valor del error: ");
    scanf("%lf", &e);

    h = (b - a) / n;

    start = clock(); // Inicio de la medición del tiempo

    for (int i = 1; i <= n; ++i) {
        yp = y + h * f(x, y);
        do {
            yc = y + 0.5 * h * (f(x, y) + f(x + h, yp));
            rel = (yp - yc) / yc;
            if (rel < 0)
                rel *= -1;
            if (rel > e)
                yp = yc;
        } while (rel > e);
        
        // Actualiza el valor de y con el valor corregido
        y = yc;
        // Incrementa x
        x = x + h;
        // Imprime los valores de x y y en cada iteración
        printf("x = %lf, y = %lf\n", x, y);
    }

    end = clock(); // Fin de la medición del tiempo

    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Tiempo de ejecucion: %lf segundos\n", cpu_time_used);
}
void ralston() {
    double x0, xf, y, h, n;
    clock_t start, end;
    double cpu_time_used;

    printf("Ingresa el valor de X: ");
    scanf("%lf", &x0); 
    
    printf("Ingresa el valor de XF: ");
    scanf("%lf", &xf); 

    printf("Ingresa el valor de Y: ");
    scanf("%lf", &y);
    
    printf("Ingresa el numero de iteraciones: ");
    scanf("%lf", &n);
    
    h = (xf - x0) / n;

    start = clock(); // Inicio de la medición del tiempo

    for (int i = 0; i < n; ++i) {
        double k1 = f(x0, y);
        double k2 = f(x0 + (3.0 / 4.0) * h, y + (3.0 / 4.0) * k1 * h);
        y = y + (h / 4.0) * (k1 + 2 * k2);
        x0 = x0 + h;
    }
    
    end = clock(); // Fin de la medición del tiempo
    
    printf("El valor de Y = %lf\n", y);
    
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Tiempo de ejecucion: %lf segundos\n", cpu_time_used);
}
void rk2() {
    double x0, xf, y, h, n;
    clock_t start, end;
    double cpu_time_used;

    printf("Ingresa el valor de X: ");
    scanf("%lf", &x0); 
    
    printf("Ingresa el valor de XF: ");
    scanf("%lf", &xf); 

    printf("Ingresa el valor de Y: ");
    scanf("%lf", &y);
    
    printf("Ingresa el numero de iteraciones: ");
    scanf("%lf", &n);
    
    h = (xf - x0) / n;

    start = clock(); // Inicio de la medición del tiempo

    for (int i = 0; i < n; ++i) {
        double k1 = f(x0, y);
        double k2 = f(x0 + h, y + k1 * h); // Corrección aquí
        y = y + (h / 2) * (k1 + k2);
        x0 = x0 + h;
    }
    
    end = clock(); // Fin de la medición del tiempo
    
    printf("El valor de Y = %lf\n", y);
    
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Tiempo de ejecucion: %lf segundos\n", cpu_time_used);
}
void heum() {
    double x0, xf, y, h, n;
    clock_t start, end;
    double cpu_time_used;

    printf("Ingresa el valor de X: ");
    scanf("%lf", &x0); 
    
    printf("Ingresa el valor de XF: ");
    scanf("%lf", &xf); 

    printf("Ingresa el valor de Y: ");
    scanf("%lf", &y);
    
    printf("Ingresa el numero de iteraciones: ");
    scanf("%lf", &n);
    
    h = (xf - x0) / n;

    start = clock(); // Inicio de la medición del tiempo

    for (int i = 0; i < n; ++i) {
        double k1 = f(x0, y);
        double k2 = f(x0 + h, y + k1 * h); // Corrección aquí
        y = y + (h / 2) * (k1 + k2);
        x0 = x0 + h;
    }
    
    end = clock(); // Fin de la medición del tiempo
    
    printf("El valor de Y = %lf\n", y);
    
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Tiempo de ejecucion: %lf segundos\n", cpu_time_used);
}
void puntomedio() {
    double x0, xf, y, h, n;
    clock_t start, end;
    double cpu_time_used;

    printf("Ingresa el valor de X: ");
    scanf("%lf", &x0); 
    
    printf("Ingresa el valor de XF: ");
    scanf("%lf", &xf); 

    printf("Ingresa el valor de Y: ");
    scanf("%lf", &y);
    
    printf("Ingresa el numero de iteraciones: ");
    scanf("%lf", &n);
    
    h = (xf - x0) / n;

    start = clock(); // Inicio de la medición del tiempo

    for (int i = 0; i < n; ++i) {
        double k1 = f(x0, y);
        double k2 = f(x0 + h/2, y + (k1 * h)/2); // Corrección aquí
        y = y + (k2 * h);
        x0 = x0 + h;
    }
    
    end = clock(); // Fin de la medición del tiempo
    
    printf("El valor de Y = %lf\n", y);
    
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Tiempo de ejecucion: %lf segundos\n", cpu_time_used);
}
void rk3() {
    double x0, xf, y, h, n;
    clock_t start, end;
    double cpu_time_used;

    printf("Ingresa el valor de X: ");
    scanf("%lf", &x0); 
    
    printf("Ingresa el valor de XF: ");
    scanf("%lf", &xf); 

    printf("Ingresa el valor de Y: ");
    scanf("%lf", &y);
    
    printf("Ingresa el numero de iteraciones: ");
    scanf("%lf", &n);
    
    h = (xf - x0) / n;

    start = clock(); // Inicio de la medición del tiempo

    for (int i = 0; i < n; ++i) {
        double k1 = f(x0, y);
        double k2 = f(x0 + h/2, y + k1 * h/2); 
        double k3 = f(x0 + h, y - k1 * h + 2 * k2 * h); 
        y = y + (h / 6) * (k1 + 4 * k2 + k3);
        x0 = x0 + h;
    }
    
    end = clock(); // Fin de la medición del tiempo
    
    printf("El valor de Y = %lf\n", y);
    
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Tiempo de ejecucion: %lf segundos\n", cpu_time_used);
}
void rk4() {
    double x0, xf, y, h, n;
    clock_t start, end;
    double cpu_time_used;

    printf("Ingresa el valor de X: ");
    scanf("%lf", &x0); 
    
    printf("Ingresa el valor de XF: ");
    scanf("%lf", &xf); 

    printf("Ingresa el valor de Y: ");
    scanf("%lf", &y);
    
    printf("Ingresa el numero de iteraciones: ");
    scanf("%lf", &n);
    
    h = (xf - x0) / n;

    start = clock(); // Inicio de la medición del tiempo

    for (int i = 0; i < n; ++i) {
        double k1 = f(x0, y);
        double k2 = f(x0 + h/2, y + k1 * h/2); 
        double k3 = f(x0 + h/2, y + k2 * h/2); 
        double k4 = f(x0 + h, y + k3 * h); 
        y = y + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
        x0 = x0 + h;
    }
    
    end = clock(); // Fin de la medición del tiempo
    
    printf("El valor de Y = %lf\n", y);
    
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Tiempo de ejecucion: %lf segundos\n", cpu_time_used);
}
void rkm5() {
    double x0, xf, y, h, n;
    clock_t start, end;
    double cpu_time_used;

    printf("Ingresa el valor de X: ");
    scanf("%lf", &x0); 
    
    printf("Ingresa el valor de XF: ");
    scanf("%lf", &xf); 

    printf("Ingresa el valor de Y: ");
    scanf("%lf", &y);
    
    printf("Ingresa el numero de iteraciones: ");
    scanf("%lf", &n);
    
    h = (xf - x0) / n;

    start = clock(); // Inicio de la medición del tiempo

    for (int i = 0; i < n; ++i) {
        double k1 = f(x0, y);
        double k2 = f(x0 + h/3, y + (h/3) * k1); 
        double k3 = f(x0 + h/3, y + (h/6) * k1 + (h/6) * k2); 
        double k4 = f(x0 + h/2, y + (h/8) * k1 + (3*h/8) * k3); 
        double k5 = f(x0 + h, y + (h/2) * k1 - (3*h/2) * k3 + 2 * h * k2);
        y = y + (h/6) * (k1 + 4*k4 + k5);
        x0 = x0 + h;
    }
    
    end = clock(); // Fin de la medición del tiempo
    
    printf("El valor de Y = %lf\n", y);
    
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Tiempo de ejecucion: %lf segundos\n", cpu_time_used);
}
void rkb6() {
    double x0, xf, y, h, n;
    clock_t start, end;
    double cpu_time_used;

    printf("Ingresa el valor de X: ");
    scanf("%lf", &x0); 
    
    printf("Ingresa el valor de XF: ");
    scanf("%lf", &xf); 

    printf("Ingresa el valor de Y: ");
    scanf("%lf", &y);
    
    printf("Ingresa el numero de iteraciones: ");
    scanf("%lf", &n);
    
    h = (xf - x0) / n;

    start = clock(); // Inicio de la medición del tiempo

    for (int i = 0; i < n; ++i) {
        double k1 = f(x0, y);
        double k2 = f(x0 + h/4, y + k1 * h / 4); 
        double k3 = f(x0 + h/4, y + (1/8) * k1 * h + (1/8) * k2 * h); 
        double k4 = f(x0 + h/2, y - (1/2) * k2 * h + k3 * h); 
        double k5 = f(x0 + (3/4) * h, y + (1/16) * k1 * h + (9/16) * k4 * h);
        double k6 = f(x0 + h, y - (3/7) * k1 * h + (2/7) * k2 * h + (12/7) * k3 * h - (12/7) * k4 * h + (8/7) * k5 * h);
        y = y + (h / 90) * (7*k1 + 32*k3 +  12*k4 + 32*k5 + 7*k6);
        x0 = x0 + h;
    }
    
    end = clock(); // Fin de la medición del tiempo
    
    printf("El valor de Y = %lf\n", y);
    
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Tiempo de ejecucion: %lf segundos\n", cpu_time_used);
}

int main(int argc, char *argv[]) {
	
	int num, op1,op2,op3,op4,opintder,opder,presder,opint, tipomod,oped;
	printf("Metodos numericos\n Proyecto Final\n");
	do{
		printf("\nMenu\n");
		printf("Elige el tema que deseas ejecutar: \n");
		printf("Temas: \n");
		printf("1.Raices de ecuaciones.\n2.Sistemas de Ecuaciones Lineales Algebraicas\n3.Derivacion e integracion numerica.\n4.Ecuaciones diferenciales ordinarias.\n5.Salir\n");
		scanf("%d",&num);
		
		switch(num){
			case 1:{
				do{
					printf("\nRaices de ecuaciones\n");
					printf("Elige el metodo a usar\n");
					printf("1.Metodo de tanteos\n2.Metodo de punto fijo\n3.Metodo de Newton-Raphson\n4.Metodo de Newton 2do orden\n5.Metodo de Biseccion\n6.Metodo de Falsa Posicion(cerrado)\n7.Metodo Secante (Hibrido)\n8.Volver\n");
					scanf("%d",&op1);
					switch(op1){
						case 1:{
							printf("Metodo de tanteos\n");
							tanteos();
							break;
						}
						case 2:{
							printf("Metodo de punto fijo\n");
							punto_fijo();
							break;
						}
						case 3:{
							printf("Metodo de Newton-Raphson\n");
							nr_1er();
							break;
						}
						case 4:{
							printf("Metodo de Newton 2do orden\n");
							newton_2do();
							break;
						}
						case 5:{
							printf("Metodo de Biseccion\n");
							biseccion();
							break;
						}
						case 6:{
							printf("Metodo de Falsa Posicion(Cerrado)\n");
							falsa_pos();
							break;
						}
						case 7:{
							printf("Metodo Secante (Hibrido)\n");
							secante();
							break;
						}
						case 8:{
							printf("Atras\n");
							break;
						}
						default:{
							printf("Opcion no valida, intenta con una opcion valida\n");
							break;
						}
					}
				}while(op1!=8);
				break;
			}
			case 2:{
				do{
					printf("Sistemas de Ecuaciones Lineales Algebraicas\n");
					printf("Elige el metodo a usar\n");
					printf("1.Algoritmo de Jacobi\n2.Algoritmo de Gauss Seidel\n3.Volver al Menu\n");
					scanf("%d",&op2);
					switch(op2){
						case 1:{
							printf("Algoritmo de Jacobi\n");
							jacobi();
							break;
						}
						case 2:{
							printf("Algoritmo de Gauss Seidel\n");
							seidell();
							break;
						}
						case 3:{
							printf("Regresando al menu\n");
							break;
						}
						default:{
							printf("Ingresa con una opcion valida\n");
							break;
						}
					}
				}while(op2!=3);
				break;
			}
			case 3:{
				do{
					printf("Derivacion e integracion numerica\n");
					printf("Elige el metodo a usar\n");
					printf("1.Derivar\n2.integrar\n3.Volver\n");
					scanf("%d",&opintder);
					switch(opintder){
						case 1:{
							do{
								printf("Derivadas\n");
								printf("1.Primera derivada\n2.Segunda Derivada\n3.Tercera derivada\n4.Cuarta derivada\n5.Volver\n");
								scanf("%d",&opder);
								switch(opder){
									case 1:{
										printf("Primera Derivada\n");
										printf("Elige la presicion de la primera derivada\n");
										printf("1.Baja Presicion\n2.Alta Presicion\n3.Volver\n");
										scanf("%d",&presder);
										do{
											switch(presder){
												case 1:{
													printf("Primera derivada Baja Presicion\n");
													pDerivdabaja();
													presder=3;
													break;
												}
												case 2:{
													printf("Primera derivada Alta Presicion\n");
													pDerivda();
													presder=3;
													break;
												}
												case 3:{
													printf("Atras\n");
													break;
												}
												default:{
													printf("Ingresa una opcion valida\n");
													break;
												}
											}
										}while(presder!=3);
										break;
									}
									case 2:{
										printf("Segunda Derivada\n");
										sDerivda();
										break;
									}
									case 3:{
										printf("Tercera Derivada\n");
										tDerivda();
										break;
									}
									case 4:{
										printf("Cuarta Derivada\n");
										cDerivada();
										break;
									}
									case 5:{
										printf("Atras\n");
										break;
									}
									default:{
										printf("Ingresa una opcion valida\n");
										break;
									}
								}
							}while(opder!=5);
							break;
						}
						case 2:{
							do{
								printf("Integrales\n");
								printf("1.Regla trapezoidal\n2.Regla de Simpson 1/3\n3.Metodo de Simpson 3/8\n4.Metodo Boole 2/45\n5.Metodo Boole 5/288\n6.Atras\n");
								scanf("%d",&opint);
								switch(opint){
									case 1:{
										do{
											printf("Regla trapezoidal\n");
											printf("Elige un modelo: \n");
											printf("1.Modelo Simple\n2. Modelo de Segmentos Multiples\n3.Salir\n");
											scanf("%d",&tipomod);
											switch(tipomod){
												case 1:{
													printf("Regla trapezoidal---Modelo Simple\n");
													trapeSimple();
													break;
												}
												case 2:{
													printf("Regla trapezoidal---Modelo de segmentos Multiples\n");
													trapemulti();
													break;
												}
												case 3:{
													printf("Atras\n");
													break;
												}
												default:{
													printf("Ingresa una opcion valida\n");
													break;
												}
											}
										}while(tipomod!=3);
									
										break;
									}
									case 2:{
										do{
											printf("Regla de Simpson 1/3\n");
											printf("Elige un modelo:\n ");
											printf("1.Modelo Simple\n2. Modelo de Segmentos Multiples\n3.Salir\n");
											scanf("%d",&tipomod);
											switch(tipomod){
												case 1:{
													printf("Regla de Simpson 1/3---Modelo Simple\n");
													simpson13Simple();
													break;
												}
												case 2:{
													printf("Regla de Simpson 1/3---Modelo de segmentos Multiples\n");
													simpson13multi();
													break;
												}
												case 3:{
													printf("Atras\n");
													break;
												}
												default:{
													printf("Ingresa una opcion valida\n");
													break;
												}
											}
										}while(tipomod!=3);
									
										break;
									}
									case 3:{
										do{
											printf("Metodo de Simpson 3/8\n");
											printf("Elige un modelo: \n");
											printf("1.Modelo Simple\n2. Modelo de Segmentos Multiples\n3.Salir\n");
											scanf("%d",&tipomod);
											switch(tipomod){
												case 1:{
													printf("Metodo de Simpson 3/8---Modelo Simple\n");
													simpson38Simple();
													break;
												}
												case 2:{
													printf("Metodo de Simpson 3/8---Modelo de segmentos Multiples\n");
													simpson38multi();
													break;
												}
												case 3:{
													printf("Atras\n");
													break;
												}
												default:{
													printf("Ingresa una opcion valida\n");
													break;
												}
											}
										}while(tipomod!=3);
									
										break;
									}
									case 4:{
										do{
											printf("Metodo Boole 2/45\n");
											printf("Elige un modelo: \n");
											printf("1.Modelo Simple\n2. Modelo de Segmentos Multiples\n3.Salir\n");
											scanf("%d",&tipomod);
											switch(tipomod){
												case 1:{
													printf("Metodo Boole 2/45---Modelo Simple\n");
													boole245Simple();
													break;
												}
												case 2:{
													printf("Metodo Boole 2/45---Modelo de segmentos Multiples\n");
													boole245multi();
													break;
												}
												case 3:{
													printf("Atras\n");
													break;
												}
												default:{
													printf("Ingresa una opcion valida\n");
													break;
												}
											}
										}while(tipomod!=3);
									
										break;
									}
									case 5:{
										do{
											printf("Metodo Boole 5/288\n");
											printf("Elige un modelo: \n");
											printf("1.Modelo Simple\n2. Modelo de Segmentos Multiples\n3.Salir\n");
											scanf("%d",&tipomod);
											switch(tipomod){
												case 1:{
													printf("Metodo Boole 5/288---Modelo Simple\n");
													boole5288Simple();
													break;
												}
												case 2:{
													printf("Metodo Boole 5/288---Modelo de segmentos Multiples\n");
													boole5288multi();
													break;
												}
												case 3:{
													printf("Atras\n");
													break;
												}
												default:{
													printf("Ingresa una opcion valida\n");
													break;
												}
											}
										}while(tipomod!=3);
									
										break;
									}
									case 6:{
										printf("Atras\n");
										break;
									}
									default:{
										printf("Ingresa una opcion valida\n");
										break;
									}
									
								}
							}while(opint!=6);
							break;
						}
						
						case 3:{
							printf("Atras\n");
							break;
						}
						default:{
							break;
						}
					}
				}while(opintder!=3);
				break;
			}
			case 4:{
				do{
					printf("Sistemas de Ecuaciones Lineales Algebraicas\n");
					printf("Elige el algoritmo a usar: \n");
					printf("1.Algoritmo de Euler \n2.Algoritmo de Euler-Gauss\n3.Algoritmo de Ralston\n4.Algoritmo de Runge-Kutta(2do orden)\n");
					printf("5.Algoritmo de Heun\n6.Algoritmo de Punto Medio\n7.Algoritmo de Runge-Kutta(3er orden)\n8.Algoritmo de Runge-Kutta(4to orden)\n");
					printf("9.Algoritmo de Runge-Kutta-Merson(5to orden)\n10.Algoritmo de Runge-Kutta-Butcher(6to orden)\n11.Atras\n");
					scanf("%d",&oped);
					switch(oped){
						case 1:{
							printf("Algoritmo de Euler\n");
							euler();
							
							break;
						}
						case 2:{
							printf("Algoritmo de Euler-Gauss\n");
							eulerGauss();
							break;
						}
						case 3:{
							printf("Algoritmo de Ralston\n");
							ralston();
							break;
						}
						case 4:{
							printf("Algoritmo de Runge-Kutta(2do orden\n)");
							rk2();
							break;
						}
						case 5:{
							printf("Algoritmo de Heun\n");
							heum();
							break;
						}
						case 6:{
							printf("Algoritmo de Punto Medio\n");
							puntomedio();
							break;
						}
						case 7:{
							printf("Algoritmo de Runge-Kutta(3er orden)\n");
							rk3();
							break;
						}
						case 8:{
							printf("Algoritmo de Runge-Kutta(4to orden)\n");
							rk4();
							break;
						}
						case 9:{
							printf("Algoritmo de Runge-Kutta-Merson(5to orden)\n");
							rkm5();
							break;
						}
						case 10:{
							printf("Algoritmo de Runge-Kutta-Butcher(6to orden)\n");
							rkb6();
							break;
						}
						case 11:{
							printf("Atras\n");
							break;
						}
						default:{
							printf("Ingresa una opcion valida\n");
							break;
						}
					}
				}while(oped!=11);
				break;
			}
			case 5:{
				printf("Fin del programa\n");
				break;
			}
			default:{
				printf("Opcion invalida, intente de nuevo\n");
				break;
			}
		}
	}while(num!=6);
	return 0;
}
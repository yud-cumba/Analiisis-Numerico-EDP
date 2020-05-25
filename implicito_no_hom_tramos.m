%ECUACION DEL CALOR POR DIFERENCIAS FINITAS IMPLICÍTAS
clear all
clc
%******************************************************
fprintf('BIENVENIDOS AL PROGRAMA PARA DESARROLLAR\nLA ECUACION DEL CALOR POR DIFERENCIAS FINITAS IMPLICÍTAS\n');
fprintf('**Ecuación no homogénea con f(x)\n que tiene una regla de correspondencia por tramos**\n');
w=input('Pulse una tecla para continuar...','s');
%******************************************************
h=input('Ingrese distancia de cada intervalo espacial:');
k=input('Ingrese distancia de cada intervalo temporal:');
a=input('Ingrese a=');
b=input('Ingrese b=');
T=input('Ingrese el tiempo final T=');
C=input('Ingrese el coeficiente C=');
m=(b-a)/h;
n=T/k;
x=linspace(a,b,m+1);
t=linspace(0,T,n+1);
r=C^2*k/(h^2);
ss=strcat('coeficiente r=',num2str(r));
disp(ss);
p=input('Ingrese punto de discontinuidad p:');
fun1=input('Ingrese la función f definida en [a,p>:','s');
fun2=input('Ingrese la fución g definida en [p,b]:','s');
f1=inline(fun1);
f2=inline(fun2);
%%%%%%%%%%%%%%%%%%%%
A=(1+2*r)*eye(m-1)+diag((-r)*ones(1,m-2),+1)+diag((-r)*ones(1,m-2),-1);
d=((p-a)/h);
%+1 si [a,p]
for i=2:d
    solfila(i-1)=f1(x(i));
end
for i=d+1:m
    solfila(i-1)=f2(x(i));
end
%solucion para t=0
h=input('Ingrese la función en la frontera a:','s');
g=input('Ingrese la función en la frontera b:','s');
alfafrontera=inline(h);
betafrontera=inline(g);
sol=[alfafrontera(t(1)),solfila,betafrontera(t(1))];
 for j=2:n+1
     H=[-r*alfafrontera(t(j)),zeros(1,m-3),-r*betafrontera(t(j))];
     fila=inv(A)*(solfila'-H');
     solfila=fila';
     sol=[sol;alfafrontera(t(j)),solfila,betafrontera(t(j))];
 end
disp('La matriz de la solución está dado por:')
u=sol';
fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
    fprintf('Para t=%.4f tenemos:\n',T);
    for i=1:m+1
        fprintf('x(%3d)=%.5f:',i-1,x(i));
        fprintf('\t %25.10f',u(i,n+1));
        fprintf('\n');
    end
 %graficat
 Tfinal=T*ones(1,m+1);
 solfinal=sol(n+1,:);
 plot3(x,Tfinal,solfinal);
 title('Solución aproximada para tiempo final por el método IMPLÍCITO');
 grid on;

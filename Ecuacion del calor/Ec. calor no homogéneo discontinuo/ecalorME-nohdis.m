%Ecuacion del Calor por Diferencias Finitas Explicitas
%MÉTODO EXPLÍCITO
clear all
clc
%******************************************************
fprintf('BIENVENIDOS AL PROGRAMA PARA DESARROLLAR\nLA ECUACION DEL CALOR POR MÉTODO EXPLÍCITO\n');
fprintf('**Ecuación no homogénea con punto de discontinuidad**\n');
w=input('Pulse una tecla para continuar...','s');
%*****************************************************
h=input('Ingrese distancia de cada intervalo espacial:');
k=input('Ingrese distancia de cada intervalo temporal:');
a=input('Ingrese el extremo a = ');
b=input('Ingrese el extremo b = ');
T=input('Ingrese el tiempo final T = ');
c=input('Ingrese el coeficiente alpha = ');
m=(b-a)/h; n=T/k;
x=linspace(a,b,m+1);
t=linspace(0,T,n+1);
r=((c^2)*k)/(h^2);
ss=strcat('coeficiente r = ',num2str(r));disp(ss); %Muestra el lambda
%Construyendo la matriz de ensamblaje
A=(1-2*r)*eye(m-1)+diag(r*ones(1,m-2),+1)+diag(r*ones(1,m-2),-1);
%Solucion inicial para t=0 en el interior de [a,b]
p=input('Ingrese el punto de discontinuidad en f: p = ');
fun1=input('Dar la funcion definida en <a,p]: f1(x) = ','s');
fun2=input('Dar la funcion definida en <p,b>: f2(x) = ','s');
f1=inline(fun1);
f2=inline(fun2);
h1un=input('Dar u(a,t) = h1(t): h1(t) = ','s');
h2un=input('Dar u(b,t) = h2(t): h2(t) = ','s');
h1=inline(h1un);
h2=inline(h2un);
xx=a+h:h:p;
xxx=p+h:h:b-h;
d=((p-a)/h)+1; %número de nodo que es p
u(1,1:d-1)=f1(xx);
u(1,d:m-1)=f2(xxx);
y=input('Pulse una tecla para continuar... ','s');
%Solucion inicial para t=0 en todo el intervalo [a,b]
sol=[h1(t(1)),u,h2(t(1))];
%El proceso iterativo en el metodo explicito
for j=2:n+1
    fila=A*u' + r*[h1(t(j-1));zeros(m-3,1);h2(t(j-1))];
    u=fila';
    sol=[sol;h1(t(j)),u,h2(t(j))];
end
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
 title('Solución aproximada para tiempo final por el método EXPLÍCITO');
 grid on;
%MÉTODO EXPLÍCITO
clear all
clc
%******************************************************
fprintf('BIENVENIDOS AL PROGRAMA PARA DESARROLLAR\nLA ECUACION DEL CALOR POR MÉTODO EXPLÍCITO\n');
fprintf('**Ecuación con punto de discontinuidad**\n');
w=input('Pulse una tecla para continuar...','s');
%******************************************************
h=input('Ingrese distancia de cada intervalo espacial:');
k=input('Ingrese distancia de cada intervalo temporal:');
a=input('Ingrese el extremo a = ');
b=input('Ingrese el extremo b = ');
T=input('Ingrese el tiempo final T = ');
c=input('Ingrese el coeficiente C = ');
m=(b-a)/h; n=T/k;
x=linspace(a,b,m+1);
tnodos=linspace(0,T,n+1);
r=((c^2)*k)/(h^2);
ss=strcat('coeficiente r = ',num2str(r));disp(ss);
%Construyendo la matriz de ensamblaje
A=(1-2*r)*eye(m-1)+diag(r*ones(1,m-2),+1)+diag(r*ones(1,m-2),-1);
%Soluci?n inicial para t=0 en el interior de [a,b]
%fun=input('Ingrese la funci?n: f(x) = ','s');
%f=inline(fun);
p=input('Ingrese el punto de discontinuidad en f: p = ');
fun1=input('Dar la funci?n definida en <a,p>: f1(x) = ','s');
fun2=input('Dar la funci?n definida en [p,b]: f2(x) = ','s');
f1=inline(fun1);
f2=inline(fun2);
alfafrontera=inline('0');
betafrontera=inline('0');
d=((p-a)/h);
%+1 [a,p]
xx=a+h:h:p-h;%p [a,p]
xxx=p:h:b-h;%p+h [a,p]
solfila(1,1:d-1)=f1(xx);
solfila(1,d:m-1)=f2(xxx);
tt=input('Pulse una tecla para continuar... ','s');
%Soluci?n inicial para t=0 en todo el intervalo [a,b]
sol=[alfafrontera(tnodos(1)),solfila,betafrontera(tnodos(1))];
%El proceso iterativo en el m?todo expl?cito
for j=2:n+1
    fila=A*solfila';
    solfila=fila';
    sol=[sol;0,solfila,0];
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
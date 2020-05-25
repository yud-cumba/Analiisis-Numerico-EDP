%Ecuacion del Calor por Diferencias Finitas Explicitas
clear all
clc
h=input('Ingrese distancia de cada intervalo espacial:');
k=input('Ingrese distancia de cada intervalo temporal:');
a=input('Ingrese a=');
b=input('Ingrese b=');

T=input('Ingrese el tiempo final T=');
c=input('Ingrese el coeficiente C=');
m=(b-a)/h;
n=T/k;
x=linspace(a,b,m+1);
t=linspace(0,T,n+1); 
r=((c^2)*k)/(h^2);
ss=strcat('coeficiente r = ',num2str(r));disp(ss); %Muestra el lambda
%Construyendo la matriz de ensamblaje
A=(1-2*r)*eye(m-1)+diag(r*ones(1,m-2),+1)+diag(r*ones(1,m-2),-1);
%Solucion inicial para t=0 en el interior de [a,b]

h1un=input('Dar u(a,t) = h1(t): h1(t) = ','s');
h2un=input('Dar u(b,t) = h2(t): h2(t) = ','s');
h1=inline(h1un);
h2=inline(h2un);

fun=input('Dar la funcion definida en t=0: f(x) = ','s');
f=inline(fun);

y=input('Pulse una tecla para continuar... ','s');
%Solucion inicial para t=0 en todo el intervalo [a,b]
for i=2:m
    u(i-1)=f(x(i));
end
sol=[h1(t(1)),u,h2(t(1))];
%El proceso iterativo en el metodo explicito
for j=2:n+1
    fila=A*u' + r*[h1(t(j-1));zeros(m-3,1);h2(t(j-1))];
    u=fila';
    sol=[sol;h1(t(j)),u,h2(t(j))];
end
u=sol';
fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
    fprintf('Para t=%.4f tenemos:\n',t(n+1));
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
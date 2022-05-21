%% Práctica 4, 13 de mayo de 2022

format long
%% 
% Problema 1

clear
a=0;
b=50;
h=0.01; %tamaño del pasos
N=(b-a)/h;
t=a:h:b;
m=1; %son iguales para cada apartado
k=1; %son iguales para cada apartado

% preparamos dos figuras para dibujar alternativamente
figure(1); clf; grid on
figure(2); clf; grid on

valores=[-2.01;-2;-1;0;1]; %van incluidos con el signo

for j=1:length(valores)
    x=zeros(2,N+1);
    c=sqrt(4*k*m)+valores(j); %te calcula c según los valores que hemos declarado
        
    % resolvemos el problema para cada valor de c
    
    f=@(t,x) [x(2); -k*x(1)/m-c*x(2)/m]; %definimos la función
    x(:,1)=[1;0]; %paso 1
    %hay que calcular los pasos del 2 al 4
    for n=1:4
        k1=f(t(n),x(:,n));
        k2=f(t(n)+h,x(:,n)+h*k1);
        x(:,n+1)=x(:,n)+h*0.5*(k1+k2);
        
    end
    
    %ahora implementamos al adam bashforth
    for n=1:N-4
        x(:,n+5)=x(:,n+4)+h*(1901/720*f(t(n+4),x(:,n+4))-2774/720*f(t(n+3),x(:,n+3))+2616/720*f(t(n+2),x(:,n+2))-1274/720*f(t(n+1),x(:,n+1))+251/720*f(t(n),x(:,n)));
    end
    
    % dibujamos en la figura 1
    figure(1)
    hold on
    plot(t,x(1,:))
    
    % dibujamos en la figura 2
    figure(2)
    hold on
    plot(x(1,:),x(2,:))

end
%azul es subamortiguado a)
%rojo el otro subamorgiuado b)
%amarillo subamortiguado c)
%magenta es el críticamente amortiguado
%verde es el super amortiguado

%queremos que el líquido no se agite, si nos fijamos en la primera gráfica
%el a) y el b) oscilan mucho
%el amarillo, el c) también tiene oscilaciones aunque menos que los otros
%entre el d) y el e) yo escogería el verde, es decir, el e) el
%superamortiguado ya que en la gráfica segunda se ve que la velocidad es
%pequeña, la verde está más cerca del eje x y a parte en la primera gráfica
%los cambios de la morada son más bruscos en comparación con la primera

%% 
% Problema 2

clear

a=0;
delta=1e-4;
b=2/delta;
h=0.1;
N=(b-a)/h;

%definimos el mallado
x=zeros(1,N+1);
t=linspace(a,b,N+1);

%la función que queremos resolver
f=@(t,x) x^2-x^3;

%valor inicial
x0=delta;
tspan=[a,b];
[d,e]=ode45(f,tspan,x0); %la 'd' es el tiempo y la 'e' es la x
figure
plot(d,e)
%en t=1 si hacemos zoom, va dando bandazos (parece un electro cardiograma), eso es porque el problema es rígido. Hay
%que intentar que la región de estabilidad sea más grande. 
%vamos a probar a resolverlo con un método de BDF de 5 pasos para
%solucionar este problema

%%

%voy a volver a poner los datos iniciales para que no se solapen con el
%apartado a
a=0;
delta=1e-4;
b=2/delta;
h=0.1;
N=(b-a)/h;

%definimos el mallado
x=zeros(1,N+1);
t=linspace(a,b,N+1);
x(1)=delta; %ponemos el dato inicial

%la función que queremos resolver
f=@(t,x) x^2-x^3;

%tenemos el primer paso x0
%hay que calcular los otros 4 con euler mejorado

for n=1:4
        k1=f(t(n),x(n));
        k2=f(t(n)+h,x(n)+h*k1);
        x(n+1)=x(n)+h*0.5*(k1+k2);        
end

%el resto de pasos con BF5

for n=1:N-4
    x(n+5)=x(n+4)+h*f(t(n+4),x(n+4)); %primera aprox con euler ya que el BDF es implícito
    for i=1:100 %punto fijo
        x_viejo=x(n+5);
        x(n+5)=300/137*x(n+4)-300/137*x(n+3)+200/137*x(n+2)-75/137*x(n+1)+12/137*x(n)+h*60/137*f(t(n+5),x(n+5));
        if abs(x_viejo-x(n+5))/abs(x(n+5))<1e-10 %el error
            break
        end
        if i==100
            fprintf('se ha llegado a 100 iteraciones') %para asegurar que no llegamos a las iteraciones máximas
        end
    end
        
end
figure
plot(t,x)

%al hacer zoom vemos que ya no va dando bandazos esto es porque el método
%BDF es más adecuado para resolver este problema rígido
% porque con este método BDF se amplia la región de
%estabilidad
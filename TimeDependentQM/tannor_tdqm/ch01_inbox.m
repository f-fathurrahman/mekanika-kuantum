% particle_in_half_box

clear;

hbar = 1;
m = 1;
L = 5;

% Number of grid points.
Nx = 2^8;
dx = 2*L/Nx;
x = [0:Nx-1]'*dx;

% Number of eigenstates included
Neigen = 30;

% Eigenstates
P = zeros(Nx,Neigen);

% Eigenstate coefficients
an = zeros(Neigen,1);

% Spectrum (Eigenvalues)
En = zeros(Neigen,1);

% Number of periods to propagate. 
Nper = 1;

% time for one period
tau = 16*m*L^2/hbar/pi;

% total propagation time.
T = tau*Nper;

Nt = 60;
dt = T/Nt;
t = [0:dt:T]';

C = zeros(Nt,1);

dw = pi/T;
wmax = Nt*dw;

% constant copmlex -i
c=-i;

%t=0;

t1 = [0.0:dt:T];
w = [0.0:dw:wmax];


% preperation of eigenstates and coeff's
for n=1:Neigen
   P(:,n) = L^-0.5 * sin(n*pi*x/2/L);
   En(n) = (n*pi*hbar)^2 / 8 / m / L^2;
   if n==2
      an(n) = 2^-0.5;
   elseif mod(n,2) == 0
      an(n) = 0;
   else
      an(n) = 4 * 2^0.5 * (-1)^((n + 1)/2) / (n+2) / (n-2) / pi;
   end
end


in = 1;
h = figure(1);
set(h,'position',[0, 50,300,250])
%set(h,'position',[0, 50,300,375])%----------------------------

psi0 = sum(P*an,2); % initial wavefunction;

% Propagation loop
for tindx = 1:Nt+1
   t = (tindx-1)*dt;  % set time
   psin = P*(an.*exp(i*En*t));
   psi = sum(psin,2);
   if tindx == Nt+1
      plot(x, abs(psi)); % plot the absolute value
      axis([0, 2*L, 0, .7]);
      f(in) = getframe(h);
      in = in + 1;   
   end

   % XXX why need this?
   if mod(tindx,1) == 0
      %f = figure('visible','off');
      %plot(x,abs(psi).^2);   
      plot(x, abs(psi));
      %plot(x,real(psi));
      axis([0, 2*L, 0, .7]);
      f(in) = getframe(h);
      in = in + 1;
   end
   %
   % trapz(x, abs(psi).^2),  % normalization check
   C(tindx) = psi0'*psi*dx;  % autocorrelation function. 
   decay = exp(-0.008*t.*t); 
   C2 = decay .* C;
end;
%movie2avi(f,'ibox.avi','compression','cinepak');

h2 = figure(2);
set(h2, 'position', [0, 50,500,400])

plot([0:Nt]*dt/tau, real(C));
title('Autocorrelation function');   %,'fontsize',20);
ylabel('C(t)','fontsize',20);
xlabel('time (in periods)','fontsize',20);
 
%spec=2*(Nt+1)*real{ifft(fftshift(C),Nt+1)}*dt-1;
for n = 1:Nt+1
   sigma1(n) = trapz(t1, C'.*exp(-c*w(n)*t1))*dt;
   sigma2(n) = trapz(t1, C2'.*exp(-c*w(n)*t1))*dt;
   spectrum(n) = 2*real(sigma1(n))-1; 
   spectrum2(n) = 2*real(sigma2(n))-1;    
end
% print -djpeg -r60 'in_2.jpeg';

h3 = figure(3);
set(h3, 'position', [0, 50,500,400])
w0 = 2*pi/tau;
w = w/w0;
plot(w, spectrum);
%xlabel('Energy (units of E_1)','fontsize',20);        
%print -djpeg -r60 'in_3.jpeg';
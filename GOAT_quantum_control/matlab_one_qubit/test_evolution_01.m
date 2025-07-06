% using struct self
clc; clear variables;

format long 

% Global variable self
global self;

% Pauli matrices
self.sx = [0 1; 1 0];
self.sy = [0 -1i; 1i 0];
self.sz = [1  0; 0  -1];
self.I  = eye(2);

% Drift Hamiltonian
self.Ho = 0.1 * self.sz;

self.num_har = 6; % ???
self.num_c = 1;
%  rng(0);

self.r  = rand(self.num_har,self.num_c);
self.w  = 0.1 ;

self.ti  = 0;
self.tf = 4;
self.steps = 1000;
self.tspan = linspace(self.ti, self.tf, self.steps);
% self.tspan = [0 8]

self.Uv  = zeros(length(self.tspan),4);
self.H = zeros(2, 2, length(self.tspan));

self.Controls = {self.sx};
self.L = {kron(eye(2),self.sx)};
self.Grad = [length(self.tspan),4 * self.num_c * self.num_har];
self.Gradients = {};
self.GradT = [];

self.U0 = eye(2);
self.UT = eye(2);

self.A = 100*rand(size(self.r))-0.5;
% self.w = 10*rand(size(self.r));
% self.w = rand(size(self.r));
% self.Bw = rand(size(self.r));

% Final state?
self.Uf = 1j* [0 1;
               1 0];

% optimization variable dumping in                        
self.iter = [];
self.state = [];
self.infidelity = [];
self.searchdir = [];
                       
self.X = (rand(self.num_har, self.num_c)-0.5);
X = self.X(:);
X  = cat(1, X, self.w);


U = eye(2);
M0 = zeros(4 + (1 * 4 * self.num_c * self.num_har + 4), 1);
M0(1:4) = U(:);

% Test evolution
opt = odeset('RelTol',1e-9,'AbsTol',1e-9,'Stats','on');
[t,M] = ode45(@(t,M) Evolution(t,M, self.A , self.w), self.tspan, M0,opt);

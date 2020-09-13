%% Eun Hwan Shin, 2020
%
% This file is the Octave implementation of Example 9.3, which solves
% the problem in Example 9.2 by elimination of the constraint. 
% 
% Mikhail, E. M. and Ackermann, F. (1976). Observations and Least Squares,
%   Thomas Y. Crowell Company, Inc., pp. 227-228.
%
% THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
% IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
% OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN
% NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
% TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%%

clear;

% number of iterations
n_iter = 6;

% Three terrestrial camera measurements from S1, S2, S3:
photo_obs = [ ...
  14.1, 6.1, 22.1;  % point A (mm)
  16.6, 7.1, 26.3]; % point B (mm)
photo_std = 0.1 ;  % mm
p = 100; % principal distance (mm)

% distance observation (m)  
dist_obs = [ ...
 10;   % S1 and S2
  8];  % S2 and S3
dist_std = 0.05; % m

Q = diag([photo_std^2*ones(6,1); dist_std^2*ones(2,1)]);

% Set approximate coordinates of point A (x1a,x2a) and point B (x2a,x2b)
% using the geometry
x2a = p * dist_obs(1) / (photo_obs(1,1) + photo_obs(1,2));
x1a = photo_obs(1,1) * x2a / p;
x2b = p * dist_obs(1) / (photo_obs(2,1) + photo_obs(2,2));
x1b = photo_obs(2,1) * x2b / p;

% Approximate coordinates in the book
% x1a = 8; x2a = 51; x1b = 7; x2b = 41;

fprintf('Approximate coordinates: \n');
fprintf('  A = (%.2f, %.2f), B = (%.2f, %.2f)\n', ...
    x1a, x2a, x1b, x2b);

% Construct linearized model: A v + B delta = f

% Partition B into B1 and B2
B1 = [-p; 0; p; 0; p; 0];
B2 = [ ...
  photo_obs(1,1), 0, 0;
  0, -p, photo_obs(2,1);
  photo_obs(1,2), 0, 0;
  0, p, photo_obs(2,2);
  photo_obs(1,3), 0, 0;
  0, p, photo_obs(2,3)];

iter = 0;
A = zeros(6,8);
A(3,7) = -p;
A(4,7) = -p;
A(5,7) = -p; A(5,8) = -p;
A(6,7) = -p; A(6,8) = -p;

while iter < n_iter,

  % increase the number of iterations
  iter = iter + 1;

  fprintf('== Iteration #%d ==\n', iter);
  
  A(1,1) = x2a;
  A(2,2) = x2b;
  A(3,3) = x2a;
  A(4,4) = x2b;
  A(5,5) = x2a;
  A(6,6) = x2b;

  f = [...
    p*x1a-photo_obs(1,1)*x2a; 
    p*x1b-photo_obs(2,1)*x2b; 
    p*dist_obs(1)-p*x1a-photo_obs(1,2)*x2a;
    p*dist_obs(1)-p*x1b-photo_obs(2,2)*x2b;
    p*dist_obs(1)+p*dist_obs(2)-p*x1a-photo_obs(1,3)*x2a;
    p*dist_obs(1)+p*dist_obs(2)-p*x1b-photo_obs(2,3)*x2b];

  % Partition constrant matrix C into C1 and C2
  C1 = 2*(x1a-x1b);
  C2 = 2*[x2a-x2b, x1b-x1a, x2b-x2a];
  g = 7.8^2 - (x1a-x1b)^2 - (x2a-x2b)^2;
  
  fprintf('Constraint: C1 = %.3f, C2 = [%.3f,%.3f,%.3f], g = %.6f\n', ...
    C1, C2(1), C2(2), C2(3), g);
 
  invC1 = 1/C1;
  B_bar = -B1 * invC1 * C2 + B2;
  f_bar = f - B1 * invC1 * g;

  QAt = Q * A';
  Qe = A * QAt;
  We = inv(Qe);
  
  N_bar = B_bar' * We * B_bar;
  t_bar = B_bar' * We * f_bar;
  
  delta2 = inv(N_bar) * t_bar;
  delta1 = invC1 * (g - C2 * delta2);

  fprintf('delta1 = %.3f, delta2 = [%.3f, %.3f, %.3f]''\n', ...
    delta1, delta2(1), delta2(2), delta2(3));

  % update parameters
  x1a = x1a + delta1;
  x2a = x2a + delta2(1);
  x1b = x1b + delta2(2);
  x2b = x2b + delta2(3);

  fprintf('Updated coordinates: \n');
  fprintf('  A = (%.3f, %.3f), B = (%.3f, %.3f)\n', ...
    x1a, x2a, x1b, x2b);

  fprintf('Distance between A and B = %.3f\n', sqrt((x1b-x1a)^2+(x2b-x2a)^2));

end	

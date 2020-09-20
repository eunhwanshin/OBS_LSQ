%% Eun Hwan Shin, 2020
%
% This file is the Octave implementation of Example 9.3, which solves
% the problem in Example 9.2 by elimination of the constraint. 
% 
% Mikhail, E. M. and Ackermann, F. (1976). Observations and Least Squares,
%   Thomas Y. Crowell Company, Inc., pp. 236-238.
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

% Problem: A v + B d = f, with Q = I
%         D1 d + D2 dp = h
% where dp denote the added parameters.

% A = eye(3);
B = [ ...
  2, -3;
 -1,  2;
  0,  1];
f = [ -1.1; 1.2; 1.0];

D1 = [ ...
  1, -1;
  2, -1];
D2 = [ 1; -2];
h = [-1; 3];

% Direct solution:
N = B' * B;
t = B' * f;
invN = inv(N);
d0 = invN * t;

P = D1 * invN * D1';
invP = inv(P);
R = D2' * invP * D2;
invR = inv(R);

dp = invR * (D2' * invP * (h - D1 * d0));
d = d0 + invN * D1' * invP * (h - D1 * d0 - D2 * dp);

% Since D1 is nonsigular and square use (9.62)
invD1 = inv(D1);
Qdd = invD1 * D2 * invR * D2' * invD1';
Qddp = -invN * D1' * invP * D2 * invR;
Qdpdp = invR;

Qxx = [ Qdd, Qddp;
  Qddp', Qdpdp];
  
% Solution by elimimation:  

D11 = D1(1,:); D12 = D2(1,:);
D21 = D1(2,:); D22 = D2(2,:);
h1 = h(1); h2 = h(2);

invD22 = inv(D22);
C = D11 - D12 * invD22 * D21;
g = h1 - D12 * invD22 * h2;

B1 = B(:,1); B2 = B(:,2);
C11 = C(1); C12 = C(2);

invC11 = inv(C11);
B_bar = -B1 * invC11 * C12 + B2;
f_bar = -B1 * invC11 * g + f;
N_bar = B_bar' * B_bar;
t_bar = B_bar' * f_bar;
invN_bar = inv(N_bar);

x2 = invN_bar * t_bar;
x1 = invC11 * (g - C12 * x2);
x3 = invD22 * (h2 - D21 * [x1; x2]);

Qx1x1 = invC11 * C12 * invN_bar * C12' * invC11;
Qx1x2 = -invC11 * C12 * invN_bar;
Qx2x2 = invN_bar;

% Verify that the two solutions are eqivalent
fprintf('d(1) - x1 = %.6f, d(2) - x2 = %.6f\n', d(1) - x1, d(2) - x2);
fprintf('dp - x3 = %.6f\n', dp - x3);
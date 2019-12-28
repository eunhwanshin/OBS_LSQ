%% Eun Hwan Shin, 2019
%
% This file is an Octave implementation of Example 8.4.
% 
% Mikhail, E. M. and Ackermann, F. (1976). Observations and Least Squares,
%   Thomas Y. Crowell Company, Inc., pp. 180-181.
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

% observations of sides of plane angles
obs = [ ...
  30, 15, 15;   % angle APB
  20,  0,  0;   % angle BPC
  50, 15, 18;   % angle APC
  30,  0,  0;   % angle CPD
  70,  0,  1];  % angle BPE
  
L = obs(:,1) + obs(:,2)/60 + obs(:,3)/3600;  

% Setup condition equations: A v + B x = f 

A = [ ...
  1, 1, -1, 0,  0;
  0, 1,  0, 1, -1];
  
B = [0; 1];  

f = -A * L;

Qe = A * A';
We = inv(Qe);

N = B' * We * B;
t = B' * We * f;

% estimate the angle DPE
Qxx = 1 / N;
x = t * Qxx;
fprintf('Angle DPE = \n');
deg2dms(x)
% residuals in arc-sec
v = A' * We * (-B * x + f) * 3600;

ref_var = v' * v;
fprintf('Reference variance = %.6f (arc-sec^2)\n', ref_var);
fprintf('var(Angle DPE) = %.6f (arc-sec^2)\n', ref_var * Qxx);

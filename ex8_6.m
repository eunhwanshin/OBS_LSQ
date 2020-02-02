%% Eun Hwan Shin, 2019
%
% This file is an Octave implementation of Example 8.6.
% 
% Mikhail, E. M. and Ackermann, F. (1976). Observations and Least Squares,
%   Thomas Y. Crowell Company, Inc., pp. 183-184.
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

% Get observations
obs_ex8_4

% observations in arc-sec
L = obs(:,1)*3600 + obs(:,2)*60 + obs(:,3);  

Q = [ ...
  2, 1, 0, 0, 0;
  1, 2, 1, 0, 0;
  0, 1, 2, 1, 0;
  0, 0, 1, 2, 1;
  0, 0, 0, 1, 2];
  
% Observations only condition: A v = f

A = [1, 1, -1, 0, 0];
f = -A * L;  

Qe = A * Q * A';
We = 1/Qe;

k = f * We;

v = Q * A' * k;

% Adjusted observations in degrees
L_hat_deg = (L+v)/3600;

B = [0, -1, 0, -1, 1];
x = B * L_hat_deg;

fprintf('Adjusted angle DPE = \n')
deg2dms(x)

% Note that with the correlation in Q, the adjusted x is different from
% what we obtained in Example 8.4 (Q = I).

Qvv = Q * A' * We * A * Q;
Qll = Q - Qvv;
Qxx = B * Qll * B';
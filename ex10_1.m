%% Eun Hwan Shin, 2020
%
% This file is the Octave implementation of Example 10.1: 
% 
% Mikhail, E. M. and Ackermann, F. (1976). Observations and Least Squares,
%   Thomas Y. Crowell Company, Inc., pp. 260-261.
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

% 6 direction measurements of a triangle has misclsure w.
% this forms an equation A v = f with Q = I, where f = w.

A = [-1, 1, -1, 1, -1, 1];

Qdd = eye(6);
Qed = A * Qdd * A';
Wed = inv(Qed);

% Thus kd = Wed*w = w/6 and the residual is
% v = A' * kd. 

Qvv = A'*Wed*A;

fprintf('6*Qvv \n')
6*Qvv

% Three intenal angle are derived from the directions by
% a = D d, where

D = [ ...
  -1,  1,  0, 0,  0, 0;
   0,  0, -1, 1,  0, 0;
   0,  0,  0, 0, -1, 1];

% Cofactor matrix for the angles
Qaa = D * Qdd * D'

% The condition equation for the anles become Aa * va = w
Aa = [1, 1, 1];

Qem = Aa * Qaa * Aa';
Wem = inv(Qem);

% thus k = w/6 and va = Qaa * Aa' * k

Qvava = Qaa * Aa' * Wem * Aa * Qaa;

fprintf('Qvava * 3/2\n')
Qvava*3/2

% Now compute Qvv from Qvava

Qll = Qdd * D' * Aa' * Wem * Aa * D * Qdd;
fprintf('Qll*6\n')
Qll*6

% This is the same as the cofactor matrix of the residual of the
% direction measuremnts. 
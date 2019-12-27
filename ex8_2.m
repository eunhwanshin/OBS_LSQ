%% Eun Hwan Shin, 2019
%
% This file is an Octave implementation of Example 8.2.
% 
% Mikhail, E. M. and Ackermann, F. (1976). Observations and Least Squares,
%   Thomas Y. Crowell Company, Inc., pp. 177-178.
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

% observations of sides of a cube (mm)
obs = [ 50.0; 50.2; 49.8]; 

% initial value of the volume
L = 49.5;
y = L^3;

iter = 0;

while iter < 5,

  iter = iter + 1;

  B = ones(3,1) * (-y^(-2/3)/3);

  f = [...
    L - obs(1);
    L - obs(2);
    L - obs(3)];

  N = B' * B;

  t = B' * f;

  delta = t/N;

  y = y + delta;
  L = y^(1/3);

  fprintf('Iteration %d: y = %.3f\n', iter, y);

end


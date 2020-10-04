%% Eun Hwan Shin, 2020
%
% This file is the Octave implementation of Example 10.2: 
% 
% Mikhail, E. M. and Ackermann, F. (1976). Observations and Least Squares,
%   Thomas Y. Crowell Company, Inc., pp. 261-263.
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

fprintf('Run Example 7.5a: adjustment of direction measurements\n')
ex7_5a 
Qvv = A_d' * We * A_d

fprintf('Run Example 7.5b: adjustment of angle measurements\n')
ex7_5b 

invQll = inv(Qll);

% assign matrix F to D
D = F;

fprintf('Rediduals in direction measurements:\n')
v_d = D' * invQll * v

DtAat = D' * A_a'
fprintf('Cofactor matrix for direction:\n')
Qvdvd = DtAat * We * DtAat'

fprintf('Note that Qvdvd derived from the angle adjustment is equivalent\n')
fprintf('to Qvv from the direction adjustement.\n')
function y = markovchain(Trans,T,s0,me,fix_seed)

% Syntax:  y=markovchain(Trans,T,s0,me,fix_seed);
% Purpose: This function generates a simulation from a Markov chain.
% INPUTS
% Trans :: transition matrix
% T     :: number of periods to be simulated
% s0    :: initial state, with one as default
% me    :: is the method (two possibilities are included, with the firstas default)

% OUTPUTS
% y     :: shock realization (indexes).
%
%Eva Carceles-Poveda 2003. Modified by Alessandro Di Nola 2018.

%Checking the mistakes from the inputs
[s1, s2]=size(Trans);

if nargin == 1
  s0=1;
  T=100;
  me=1;
  fix_seed=0;
end
if nargin == 2
   s0=1;
   me=1;
   fix_seed=0;
end
if nargin == 3
   me=1;
   fix_seed=0;
end
if nargin == 4
   fix_seed=0;
end

%check that Trans and s0 are well defined
if s1 ~= s2
  disp('Transition matrix must be square');
  return
end

for k=1:s1
  if sum(Trans(k,:)) ~= 1
    disp(['Row ',num2str(k),' does not sum to one']);
    disp(['Normalizing row ',num2str(k),'']);
    Trans(k,:)=Trans(k,:)/sum(Trans(k,:));
  end
end

if s0 < 1 ||s0 > s1
  disp(['Initial state ',num2str(s0),' is out of range']);
  disp(['Initial state defaulting to 1']);
  s0=1;
end

%Creating the shock realizations
if fix_seed == 1
    rng('default');
end
X=rand(T,1);
y(1,1)=s0;
if me==1%use method 1
   for i=2:T
      for j=1:s1   
         if X(i-1)<sum(Trans(y(i-1),1:j)) 
            break
      		end
        j=j+1;
      end
      y(i,1)=j;
   end
elseif me==2
  	s=zeros(s1,1);
	s(s0)=1;
	cum=Trans*triu(ones(size(Trans)));
	for k=1:length(X)
  		state(:,k)=s;
  		ppi=[0 s'*cum];
  		s=((X(k)<=ppi(2:s1+1)).*(X(k)>ppi(1:s1)))';
	end
   y=[1:s1]*state;
end

 
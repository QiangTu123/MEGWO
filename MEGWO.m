%！！！！！！！！！！！！！！！！！！！！！！！！！！！！------------------%
%
% Multi-strategy Ensemble Grey Wolf Optimizer 
%       By: Mr.Tu
%       Email:tuqiang3@mail2.sysu.edu.cn
%
% Developed in MATLAB R2014a
% The original code of GWO is availble on 
%     Homepage: http://www.alimirjalili.com
%
%！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！%

function [Alpha_score,Alpha_pos,Cbest]=MEGWO(SearchAgents_no,Max_iter,lb,ub,dim,f)

%----------Initialize the input parameter-----------------%
rand('state',sum(100*clock));
d=dim; 
N=SearchAgents_no;
M=Max_iter;
%---------------------------------------------------------%

% initialize alpha, beta, and delta_pos
Alpha_pos=zeros(1,d);
Alpha_score=inf; %change this to -inf for maximization problems
Beta_pos=zeros(1,d);
Beta_score=inf; %change this to -inf for maximization problems
Delta_pos=zeros(1,d);
Delta_score=inf; %change this to -inf for maximization problems
%Initialize the positions of search agents
aa=repmat(lb,N,1);
bb=repmat(ub,N,1);
Positions=aa+(bb-aa).*rand(N,d);
X=Positions;
X(:,d+1)=feval(f, X(:,1:d));
fit=X(:,d+1);
for i=1:N
    if fit(i)<Alpha_score
        Alpha_score=fit(i); % Update alpha
        Alpha_pos=Positions(i,:);
    end
    if fit(i)>Alpha_score && fit(i)<Beta_score
        Beta_score=fit(i); % Update beta
        Beta_pos=Positions(i,:);
    end
    if fit(i)>Alpha_score && fit(i)>Beta_score && fit(i)<Delta_score
        Delta_score=fit(i); % Update delta
        Delta_pos=Positions(i,:);
    end
end
a1=2.0;
l=0;% Loop counter
Cbest=zeros(1,M);
% Main loop
while l<M
    pa=1-(1-0.6)*l/M;
    freq=0.4-(0.4-0)*l/M;
    a=2*(1-l/M); % a decreases linearly fron 2 to 0
    % Update the Position of search agents including omegas    
    Positions=X(:,1:d);    
    for i=1:N
        if rand<0.5
            for j=1:d
                cnum=ceil(d*rand);
                if rand<=0.8
                    Positions(i,j)=Alpha_pos(cnum);
                else
                    rnum=ceil(N*rand);while rnum==i,rnum=ceil(N*rand);end;
                    rnum1=ceil(N*rand);while rnum1==i || rnum1==rnum,rnum1=ceil(N*rand);end;
                    Positions(i,j)=Alpha_pos(j)+2*a*rand*(Positions(rnum,j)-Positions(rnum1,j));
                end
            end
        else
            if rand<=pa
            k=ceil(rand*d);
            else
            k=1:d;   
            end
            for j=k            
                X1=Alpha_pos(j)-a*(2*rand-1)*abs(a1*rand*Alpha_pos(j)-Positions(i,j));
                X2=Beta_pos(j)-a*(2*rand-1)*abs(a1*rand*Beta_pos(j)-Positions(i,j));
                X3=Delta_pos(j)-a*(2*rand-1)*abs(a1*rand*Delta_pos(j)-Positions(i,j));
                Positions(i,j)=(X1+X2+X3)/3;
            end
        end        
    end
    Positions=ControlBoundD(Positions,aa,bb);
    fit=feval(f, Positions);
    for i=1:N        
        fitness=fit(i);
        if X(i,d+1)>fitness
            X(i,d+1)=fitness;X(i,1:d)=Positions(i,:);
        end        
        % Update Alpha, Beta, and Delta
        if fitness<Alpha_score
            Alpha_score=fitness; % Update alpha
            Alpha_pos=Positions(i,:);
        end
        if fitness>Alpha_score && fitness<Beta_score
            Beta_score=fitness; % Update beta
            Beta_pos=Positions(i,:);
        end
        if fitness>Alpha_score && fitness>Beta_score && fitness<Delta_score
            Delta_score=fitness; % Update delta
            Delta_pos=Positions(i,:);
        end        
    end       
    nest=X(:,1:d);
    K=rand(N,d)>freq;
    stepsize=normrnd(0.1,0.5)*(nest(randperm(N),:)-nest(randperm(N),:));
    Positions=nest+stepsize.*K;
    Positions=ControlBoundD(Positions,aa,bb);
    fit=feval(f, Positions);
    for i=1:N
        fitness=fit(i);
        if X(i,d+1)>fitness
            X(i,d+1)=fitness;X(i,1:d)=Positions(i,:);
        end        
        % Update Alpha, Beta, and Delta
        if fitness<Alpha_score
            Alpha_score=fitness; % Update alpha
            Alpha_pos=Positions(i,:);
        end
        if fitness>Alpha_score && fitness<Beta_score
            Beta_score=fitness; % Update beta
            Beta_pos=Positions(i,:);
        end
        if fitness>Alpha_score && fitness>Beta_score && fitness<Delta_score
            Delta_score=fitness; % Update delta
            Delta_pos=Positions(i,:);
        end
    end  
    l=l+1; 
    Cbest(l)=Alpha_score;       
end
function V=ControlBoundD(V,L,U)
c=V<L;
if sum(c(:))
    if rand<0.2
       V(c)=L(c); 
    else 
    V(c)=L(c)+(U(c)-L(c)).*rand(sum(c(:)),1);
    end
end 
c=V>U;
if sum(c(:))
    if rand<0.2
       V(c)=U(c); 
    else 
    V(c)=L(c)+(U(c)-L(c)).*rand(sum(c(:)),1);
    end
end 
function [q,X,S,NewInd]=AddOneElement(Omega,Q,X,S,Y,opt)
%AddOneElement adds one entry to the co-support, updates the signal 
%estimate and the accumulated orthogonal set.
%  [q,X,S,NewInd]=AddOneElement(Omega,Q,X,S,Y,opt)
%  ========================================================================
%  Input:
%  Omega - analysis dictionary.
%  Q     - matrix of the orthogonal set accumulated so far (in its rows).
%  X     - the current estimate for the signal.
%  S     - the current estimate for the co-support.
%  Y     - the observed (possibly noisy) signal.
%  opt   - the type of optimality criterion used for adding elements to the 
%          co-support: 0 - THR, 1 - BG, 2 - OBG. 
%  ========================================================================
%  Output:
%  q      - the new row vector added to the orthogonal set.
%  X      - the new estimate for the signal.
%  S      - the new estimate for the co-support.
%  NewInd - the new element added to the co-support.
%  ========================================================================
%  Tomer Peleg
%  Department of Electrical Engineering
%  Technion, Haifa 32000 Israel
%  tomerfa@tx.technion.ac.il
%
%  October 2012
%  ========================================================================
d=size(Omega,2);
inds=find(S==0);
q=[];
NewInd=[];
if ~isempty(inds)
    Q_temp=Omega(inds,:);
    if opt==2
        if ~isempty(Q)
            % modified gram-schmidt for all possible directions
            for j=1:size(Q,1)
                Q_temp=Q_temp-(Q_temp*Q(j,:)')*Q(j,:);
            end
            Q_temp=Q_temp./repmat(sqrt(diag(Q_temp*Q_temp')),1,d);
        end
    end
    % project on each of the possible directions
    if opt==0
        Proj=(Q_temp*Y).^2;
    else
        Proj=(Q_temp*X).^2;
    end
    % choose the direction with the smallest energy reduction
    [dummy,pos] = min(Proj);
    q=Q_temp(pos,:);
    if opt<2
        if ~isempty(Q)        
            for j=1:size(Q,1)
                q=q-(q*Q(j,:)')*Q(j,:);
            end
            q=q/sqrt(q*q');
        end
    end
    if nargout>1
        X=X-q'*(q*X);
    end
    if nargout>2
        S=abs(Omega*X)<1e-4;
    end
    if nargout>3
        NewInd=inds(pos);
    end
end
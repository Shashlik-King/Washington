function [var_q]=interp_var(depth,var,depth_q, Flag)

% depth_q=4.3;
% var=[Gmax,Gmax,Gmax];

% if strcmpi(Flag,'top')
%     I=1;
% elseif strcmpi(Flag,'')
%     I=1;
% elseif strcmpi(Flag,'bottom')
%     I=2;
% end

depth=round(depth,5);
depth_q=round(depth_q,5);


depth_q=sort(abs(depth_q));

% if length(unique(depth_q))<length(depth_q)
%     if ~strcmpi(Flag,'')
%     error('When interpolating at repeated values of depth. Flag should be empty !!!');
%     end
% end

var_q=zeros(length(depth_q),size(var,2));

for ii=1:length(depth_q)
    
    %disp(num2str(ii))
    
    x = depth_q(ii);
    
    idx0=find(depth==x);
    
    if ~isempty(idx0)
        
        if length(idx0)>1
            
            try
                x_previous = depth_q(ii-1);
            catch
                x_previous=[];
            end
            
            
            if x==x_previous % if querried depth was repeated before, then use the value at the bottom of the layer
                idx0=idx0(2);
            else
                idx0=idx0(1);
                
            end
            
        end
        
        %v=mean(var(idx0,:),1);
        v=var(idx0,:);
        %disp('Interplated point lies at an interface. An average was used!')
        
    else
        
        idx1=find(depth>x,1);
        
        if isempty(idx1)
            warning('Queried depth does not exist in the profile. Last value was used to extrapolate');
            idx1 = length(depth);
        end
        
        
        x1=depth(idx1);
        x2=depth(idx1-1);
        
        v1=var(idx1,:);
        v2=var(idx1-1,:);
        
        v = (v1-v2)/(x1-x2)*(x-x1)+v1;
        
        
    end
    
    var_q(ii,:)=v;
    
end




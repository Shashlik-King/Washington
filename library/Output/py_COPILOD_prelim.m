for i=1:length(element.level)-1
    output.excel(i*2-1,1)=element.level(i,1);
    output.excel(i*2-0,1)=element.level(i,1);
    output.excel(i*2-1,2:18)=p.top(i,:);
    output.excel(i*2-0,2:18)=y.top(i,:);
end
output.excel(i*2+1,1)=element.level(i+1,1);
output.excel(i*2+2,1)=element.level(i+1,1);
output.excel(i*2+1,2:18)=p.bottom(i,:);
output.excel(i*2+2,2:18)=y.bottom(i,:);
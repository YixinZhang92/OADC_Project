fig=figure; 
hax=axes; 
x=0:0.1:10; 
hold on 
plot(x,sin(x)) 
SP=1; %your point goes here 
line([SP SP],get(hax,'YLim'),'Color',[1 0 0])
line([2*SP 2*SP],get(hax,'YLim'),'Color',[1 1 0])
line([3*SP 3*SP],get(hax,'YLim'),'Color',[0 1 0])
line([4*SP 4*SP],get(hax,'YLim'),'Color',[0 0 1])
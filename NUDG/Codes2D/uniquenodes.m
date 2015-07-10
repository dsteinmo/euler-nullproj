%Given an Mx2 input array of nodes, return unique nodes with mapping indices i and j such that nodes(i) = unodes and unodes(j) = nodes. Depends on 'consolidator' matlab script

function [unodes,i,j] = uniquenodes(x,y)

   nodetol = 1e-7;
   dx = abs(diff(x));
   dy = abs(diff(y)); 
   dx = dx(dx > nodetol);
   dy = dy(dy > nodetol);

   min_d = min([min(dx) min(dy)]);
   min_d = min_d*(1-nodetol); %need to fudge this, or else consolidator can bork.
                         %this 0.99 seems fragile. me no like.

   min_d
   nodes = [x(:) y(:)];
   [unodes,blah,j] = consolidator(nodes,[],[],.5*min_d);

   i=zeros(size(unodes),1);

   for ind=1:length(unodes)
      point = [unodes(ind,1) unodes(ind,2)];
      pointrm = repmat(point,length(nodes),1);
      dists = sqrt((nodes(:,1)-pointrm(:,1)).^2 + (nodes(:,2)-pointrm(:,2)).^2);
      foundInd = find(dists < min_d);
      i(ind) = foundInd(1);
   end

end % uniquenodes

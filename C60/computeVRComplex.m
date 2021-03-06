function [kSkeletonOfVRComplex,simplexDimension]=computeVRComplex(G,k)
clear kSkeletonOfVRComplex
nPoints=size(G,1);
simplexDimension=k;
%1) Intialize Complex with 0-Simplex
kSkeletonOfVRComplex{1}=zeros(nPoints,1);
for v=1:size(G,1)
   kSkeletonOfVRComplex{1}(v)=v; 
end

%2) Compute k-complexes using the previous simplex set (Induction Method:
%Zomorodian, A. (2010). Fast construction of the Vietoris-Rips complex. Computers & Graphics, 34(3), 263?271. doi:10.1016/j.cag.2010.03.007)
for i=1:k
   kSimplexCount=1;
   kSimplex=[];
   nSimplicies=size(kSkeletonOfVRComplex{i},1);
   for t=1:nSimplicies
         nElementsInPreviousSimplex=size(kSkeletonOfVRComplex{i},2);
         %Intersect the lower neighbourhoods for each u in simplex T.
%          intersectionOfNeighbourhoods=lowerNeighbours(G,kSkeletonOfVRComplex{i}(t,1))';
%          for u=2:nElementsInPreviousSimplex
%              intersectionOfNeighbourhoods=intersectionOfNeighbourhoods(ismembc(intersectionOfNeighbourhoods,lowerNeighbours(G,kSkeletonOfVRComplex{i}(t,u))));
%          end
         % check if all element in the row equal to 1, if so record it
         intersectionOfNeighbourhoods=kSkeletonOfVRComplex{1}(all(G(:,kSkeletonOfVRComplex{i}(t,:))==1,2));
         
         %Create k-simplicies for each element of
         %InterserctionOfNeighbourhoods
        if ~isempty(intersectionOfNeighbourhoods) 
            for v=1:length(intersectionOfNeighbourhoods)
                
              %Simplex is the elements in the current order complex unioned with the
              %intersection of neighbours. (Ignore empty entries)
                kSimplex(kSimplexCount ,:)=sort([kSkeletonOfVRComplex{i}(t,:) intersectionOfNeighbourhoods(v)]) ;
                kSimplexCount=kSimplexCount+1;
            end
        end
   end
        if isempty(kSimplex)
            simplexDimension=i;
            for m=i+1:k+1
                kSkeletonOfVRComplex{m}=[];
            end
            return
        end
            

    %Add the set of k complexes as another cell array in the skeleton.
    kSkeletonOfVRComplex{i+1}=unique(kSimplex,'rows');

end

end


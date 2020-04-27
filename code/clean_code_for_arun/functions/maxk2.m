function [B BIndex RestVector]= maxk2(A, k)
B = 0;
RestVector = A;
sumIndex = 1;
for i=1:k
  MaxA = max(A);
  I = A == MaxA;
  sumI = sum(I); %To find number of Max elements (repeated) 
  B(sumIndex: sumIndex+sumI-1) = MaxA; % to same max elements in B
  BIndex(sumIndex: sumIndex+sumI-1) = find(A == MaxA); 
  sumIndex = sumIndex + sumI; 
  A(I) = min(A); % exchange the max elements by a smallest value  
end
RestVector(BIndex) = [];  % remove largest values

end


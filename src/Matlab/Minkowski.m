function M=minksum(P,Q)

%computes the minkowski sum of two convex polygons: P and Q. The polygons
%are represented by their vertices and are ordered counter clockwise such
%that the first vertex will be the one who has the smallest Y coordinate
%(and smallest X coordinate in case of a tie).

m=size(P,1);
n=size(Q,1); %number of verteces

i=1; j=1;
M=[]; %M will contain the Minkowski sum

% it is important that the P polygon will end with the largest angle and 
% not the Q polygon. That means that if we look at the last vertex of P and
% the last vertex of Q and compare the angle that each creates with the 
% first vertex of the corresponding polygon. These angles are denoted as
% angP and angQ. If angP<angQ then we need to switch between them.

tol=1E-10;
angP=xangle(P(m,1),P(m,2),P(1,1),P(1,2));
angQ=xangle(Q(n,1),Q(n,2),Q(1,1),Q(1,2));

if angP>angQ
    while (i<m+1) | (j<n+1)
        current_i=min(mod(i,m+1)+1,i);
        current_j=min(mod(j,n+1)+1,j);
        next_i=min(mod(i+1,m+1)+1,i+1);
        next_j=min(mod(j+1,n+1)+1,j+1);
        M=[M; P(current_i,:)+Q(current_j,:)];
        angP=xangle(P(current_i,1),P(current_i,2),P(next_i,1),P(next_i,2));
        angQ=xangle(Q(current_j,1),Q(current_j,2),Q(next_j,1),Q(next_j,2));
        %[i j angP angQ]
        iinc=0;
        jinc=0;
        if (angP<=angQ+tol) | (j==n+1)
            iinc=1;
        end
        if (angP>=angQ-tol) | (i==m+1)
            jinc=1;
        end
        i=i+iinc;
        j=j+jinc;
        i=min(i,m+1);
        j=min(j,n+1);
    end
    %    PP=[P; P(1:2,:)];
    %    QQ=[Q; Q(1:2,:)]; %we add the first vertex at the end of the list
else %lets switch between P and Q
    while (i<m+1) | (j<n+1)
        current_i=min(mod(i,m+1)+1,i);
        current_j=min(mod(j,n+1)+1,j);
        next_i=min(mod(i+1,m+1)+1,i+1);
        next_j=min(mod(j+1,n+1)+1,j+1);
        M=[M; Q(current_j,:)+P(current_i,:)];
        angP=xangle(P(current_i,1),P(current_i,2),P(next_i,1),P(next_i,2));
        angQ=xangle(Q(current_j,1),Q(current_j,2),Q(next_j,1),Q(next_j,2));
        %[i j angP angQ]
        iinc=0;
        jinc=0;
        if (angQ<=angP) | (i==m+1)
            jinc=1;
        end
        if (angQ>=angP) | (j==n+1)
            iinc=1;
        end
        i=i+iinc;
        j=j+jinc;
        i=min(i,m+1);
        j=min(j,n+1);
    end

%    PP=[Q; Q(1:2,:)];
%    QQ=[P; P(1:2,:)]; %we add the first vertex at the end of the list
end
% while (i<m+1) | (j<n+1)
%     current_i=min(mod(i,m+1)+1,i);
%     current_j=min(mod(j,n+1)+1,j);
%     next_i=min(mod(i+1,m+1)+1,i+1);
%     next_j=min(mod(j+1,n+1)+1,j+1);
%     M=[M; P(current_i,:)+Q(current_j,:)];
%     angP=xangle(P(current_i,1),P(current_i,2),P(next_i,1),P(next_i,2));
%     angQ=xangle(Q(current_j,1),Q(current_j,2),Q(next_j,1),Q(next_j,2));
%     %[i j angP angQ]
%     if (angP<=angQ) | (j==n+1)
%         i=i+1;
%     end
%     if (angP>=angQ) | (i==m+1)
%             j=j+1;
%     end
%     i=min(i,m+1);
%     j=min(j,n+1);
% end
%[i j]
%M=[M; PP(m,:)+QQ(n,:)];
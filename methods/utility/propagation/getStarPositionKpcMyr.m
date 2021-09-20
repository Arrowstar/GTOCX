function [rVectKpc, vVectKpcMyr] = getStarPositionKpcMyr(ids, tMyr, starData)
%getStarPosition Returns the position and velocity of given stars at given
%times.
%   INPUTS:
%       ids - a column vector of star ID numbers as given in stars.txt.
%       tMyr - a column vector of times to compute the star
%       positions/velocities, given in Myr
%       starData - the matrix of data from the stars.txt file
%   OUTPUTS:
%       rVectKpc - a 3 row, N column matrix of star positions.  Each column
%       is a position vector in kpc.
%       vVectKpcMyr - a 3 row, N column matrix of star velocities.  Each
%       column is a star velocity vector in kpc/Myr

%     rowNum = ids+1;
    ids = ids(:);
    tMyr = tMyr(:);

    [~,Locb] = ismember(ids,starData(:,1));
    rowData = starData(Locb,:);

    Rkpc = rowData(:,2);
    i = deg2rad(rowData(:,3));
    Om = deg2rad(rowData(:,4));
    phi = deg2rad(rowData(:,5));
    
    VcKpcMyr = computeVc(Rkpc);                                            %(1)
    
    n = VcKpcMyr./Rkpc; %1/Myr                                              (6)
    t = tMyr; %Myr
    
    x = Rkpc .* (cos(n.*t+phi).*cos(Om) - sin(n.*t+phi).*cos(i).*sin(Om));  %(7)
    y = Rkpc .* (cos(n.*t+phi).*sin(Om) + sin(n.*t+phi).*cos(i).*cos(Om));  %(8)
    z = Rkpc .* (sin(n.*t+phi).*sin(i));                                    %(9)
    
    vx = VcKpcMyr .* (-sin(n.*t + phi).*cos(Om) - cos(n.*t + phi).*cos(i).*sin(Om)); %(10)
    vy = VcKpcMyr .* (-sin(n.*t + phi).*sin(Om) + cos(n.*t + phi).*cos(i).*cos(Om)); %(11)
    vz = VcKpcMyr .* (cos(n.*t + phi).*sin(i));                                      %(12)
    
    rVectKpc = [x';y';z']; %kpc
    vVectKpcMyr = [vx';vy';vz'];
end


function r = buildCoordRotationDataStraight
%like buildCoordRotationData but gives data for only
%elements that are straight
    global straight
    
    r=buildCoordRotationData;
    r.A11 = r.A11(:,straight);
    r.A12 = r.A12(:,straight);
    r.A21 = r.A21(:,straight);
    r.A22 = r.A22(:,straight);
    
    r.Ainv11 = r.Ainv11(:,straight);
    r.Ainv12 = r.Ainv12(:,straight);
    r.Ainv21 = r.Ainv21(:,straight);
    r.Ainv22 = r.Ainv22(:,straight);
    
    r.x0 = r.x0(:,straight);
    r.y0 = r.y0(:,straight);

end
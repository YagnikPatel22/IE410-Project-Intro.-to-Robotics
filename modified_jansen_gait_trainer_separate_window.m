%% modified_jansen_gait_trainer.m
% Kinematic reproduction of the single-DOF modified Jansen gait trainer.
%
% The model follows the 12-length, eight-bar topology used in:
% Shin, Deshpande, Sulzer, JMR 2018, and the later validationGrid paper.
% Units are cm, seconds, and radians unless otherwise noted.
%
% Run this file. It generates:
%   1. A paper-style Fig. 4 output: mechanism, endpoint path, x/y gait curves.
%   2. A paper-style 3x3 validationGrid output with RMSE values.
%   3. A running animation of the one-DOF mechanism.
%   4. Optional diagnostic plots and L1/L4/L8 span analysis.

clear; close all; clc;

if exist("fsolve", "file") == 2
    fprintf("Solver: MATLAB fsolve from Optimization Toolbox.\n");
else
    fprintf("Solver: included damped Newton nonlinearSolver because fsolve is unavailable.\n");
end

%% 1. Parameters
% Nominal dimensions from the validationGrid paper appendix, with the paper's
% example adjustable values substituted for L1, L4, and L8.
% Appendix nominal set: [11,45,36,34,48.5,41.5,60.5,41.5,42,43,26.5,54.5].
linkLen = struct();
linkLen.names = ["L1","L2","L3","L4","L5","L6","L7","L8","L9","L10","L11","L12"];
linkLen.val = [11.29, 45.0, 36.0, 32.93, 48.5, 41.5, 60.5, 41.78, 42.0, 43.0, 26.5, 54.5].';

cfg = struct();
cfg.N = 361;                    % one full crankAngles turn
cfg.omega = 2*pi/2.0;           % rad/s; 2 s per gait cycle
cfg.phaseOffset = 0.0;          % radians; changes starting point only
cfg.gaitFrameRotationDeg = -11.4; % rotates mechanism frame to gait/world frame
cfg.showPaperFigure4 = true;
cfg.showNinePatternValidation = true;
cfg.runAnimation = true;
cfg.animationFrameStep = 3;      % higher value makes animation faster
cfg.showDiagnosticPlots = false; % set true for separate trajectory/velocity/residual figures
cfg.showMechanismFrames = false;
cfg.NAdvanced = 91;             % lower resolution for sensitivityData sweeps
cfg.runSensitivity = false;      % set true for adjustable-link sensitivityData plots
cfg.runLeastSquaresDemo = false; % set true for span-to-link least-squares demo

fprintf("Modified Jansen gait trainer simulation\n");
fprintf("Using adjustable values: L1 = %.2f cm, L4 = %.2f cm, L8 = %.2f cm\n", ...
    linkLen.val(1), linkLen.val(4), linkLen.val(8));

%% 2. Main one-cycle simulation
trajData = simulateJansenCycle(linkLen.val, cfg.N, cfg.omega, cfg.phaseOffset, [], ...
    cfg.gaitFrameRotationDeg*pi/180);

spanX = spanOf(trajData.ptE(1,:));
spanY = spanOf(trajData.ptE(2,:));
areaPE = abs(polyarea(trajData.ptE(1,:), trajData.ptE(2,:)));
fprintf("End-effector x-span = %.2f cm\n", spanX);
fprintf("End-effector y-span = %.2f cm\n", spanY);
fprintf("Closed-loop area     = %.2f cm^2\n", areaPE);
fprintf("Max loop residual    = %.3e cm\n\n", max(trajData.resnorm));

%% 3. Paper-style reference curve for visual comparison
% The actual human marker data used in the papers is not inside the PDFs.
% This pchip curve recreates the reported gait envelope:
% desired span [xspan, yspan] = [50.02, 12.81] cm.
targetSpan = [50.02; 12.81];
refCurve = makeReferenceGait(trajData.gaitPercent, targetSpan);
refAligned = alignCurveToSimulation(refCurve, trajData.ptE);

%% 4. Paper-style outputs
if cfg.showPaperFigure4
    plotPaperStyleFigure4(trajData, refAligned, linkLen);
end

if cfg.showNinePatternValidation
    validationGrid = plotNinePatternValidation(linkLen.val, cfg);
end

if cfg.runAnimation
    animateJansenMechanism(trajData, linkLen, cfg.animationFrameStep);
end

%% 5. Optional diagnostic plots
if cfg.showDiagnosticPlots
    plotEndEffector(trajData, refAligned, targetSpan, linkLen);
    plotGaitCurves(trajData, refAligned);
    plotVelocityCurves(trajData);
    plotConstraintResiduals(trajData);
end

if cfg.showMechanismFrames
    plotMechanismSnapshots(trajData, linkLen);
end

%% 6. Optional advanced plots: sensitivityData and least-squares link-span mapping
if cfg.runSensitivity
    sensitivityData = runAdjustableLinkSensitivity(linkLen.val, cfg);
    plotSensitivity(sensitivityData);
end

if cfg.runLeastSquaresDemo
    lsMap = runLeastSquaresSpanMap(linkLen.val, cfg, targetSpan);
    plotLeastSquaresDemo(lsMap, targetSpan);
end

%% 7. Data exported to the MATLAB workspace
result = struct();
result.lengths_cm = linkLen;
result.trajData = trajData;
result.reference = refAligned;
result.targetSpan_cm = targetSpan;
if exist("validationGrid", "var")
    result.validationGrid = validationGrid;
end
if exist("sensitivityData", "var")
    result.sensitivityData = sensitivityData;
end
if exist("lsMap", "var")
    result.leastSquaresDemo = lsMap;
end

disp("Done. The struct named 'result' contains trajectories, angles, spans, and residuals.");

%% Local functions

function cycleState = simulateJansenCycle(linkLen, N, omega, phaseOffset, seedAngles, frameRot)
%SIMULATEJANSENCYCLE Solve all closed loops for one crankAngles revolution.
%
% Link vector directions used in this code:
%   L1:  pt0 -> pt1       L2:  pt1 -> pt2       L3:  pt3 -> pt2
%   L4:  pt0 -> pt3       L5:  pt2 -> pt4       L6:  pt3 -> pt4
%   L7:  pt1 -> pt5       L8:  pt3 -> pt5       L9:  pt5 -> pt6
%   L10: pt4 -> pt6       L11: pt5 -> ptE       L12: ptE -> pt6
%
% Unknown vector:
%   unknownAngles = [theta2 theta3 theta5 theta6 theta7 theta8 theta9 theta10 theta11 theta12]
% theta1 is the prescribed input crankAngles angle; theta4 = 0 is the ground angle.

    if nargin < 6
        frameRot = 0;
    end

    crankAngles = linspace(0, 2*pi, N);
    timeAxis = crankAngles / omega;
    gaitPercent = 100 * crankAngles / (2*pi);

    nonlinearSolver = makeLoopSolver();

    angleMat = nan(12, N);
    qHist = nan(10, N);
    exitflag = nan(1, N);
    resnorm = nan(1, N);

    pt0 = nan(2, N); pt1 = pt0; pt2 = pt0; pt3 = pt0; pt4 = pt0;
    pt5 = pt0; pt6 = pt0; ptE = pt0;

    unknownAngles = seedAngles;
    for k = 1:N
        theta1 = crankAngles(k) + phaseOffset;

        if isempty(unknownAngles) || any(~isfinite(unknownAngles))
            unknownAngles = geometricInitialGuess(linkLen, theta1);
        end

        loopFun = @(qq) loopEquations(qq, theta1, linkLen);
        loopJac = @(qq) loopJacobian(qq, linkLen);
        [qCandidate, fval, flag] = solveLoopAngles(loopFun, loopJac, unknownAngles, nonlinearSolver);
        if flag <= 0 || norm(fval) > 1e-6
            qGeom = geometricInitialGuess(linkLen, theta1);
            [qCandidate2, fval2, flag2] = solveLoopAngles(loopFun, loopJac, qGeom, nonlinearSolver);
            if norm(fval2) < norm(fval)
                qCandidate = qCandidate2;
                fval = fval2;
                flag = flag2;
            end
        end

        unknownAngles = unwrapNear(qCandidate(:), unknownAngles(:));
        qHist(:,k) = unknownAngles;
        exitflag(k) = flag;
        resnorm(k) = norm(fval);

        linkAngles = nan(12,1);
        linkAngles(1) = theta1;
        linkAngles(2) = unknownAngles(1);
        linkAngles(3) = unknownAngles(2);
        linkAngles(4) = 0;
        linkAngles(5) = unknownAngles(3);
        linkAngles(6) = unknownAngles(4);
        linkAngles(7) = unknownAngles(5);
        linkAngles(8) = unknownAngles(6);
        linkAngles(9) = unknownAngles(7);
        linkAngles(10) = unknownAngles(8);
        linkAngles(11) = unknownAngles(9);
        linkAngles(12) = unknownAngles(10);
        angleMat(:,k) = linkAngles;

        pointSet = positionsFromAngles(linkLen, linkAngles);
        pt0(:,k) = pointSet.pt0; pt1(:,k) = pointSet.pt1; pt2(:,k) = pointSet.pt2;
        pt3(:,k) = pointSet.pt3; pt4(:,k) = pointSet.pt4; pt5(:,k) = pointSet.pt5;
        pt6(:,k) = pointSet.pt6; ptE(:,k) = pointSet.ptE;
    end

    raw = struct("pt0", pt0, "pt1", pt1, "pt2", pt2, "pt3", pt3, "pt4", pt4, ...
        "pt5", pt5, "pt6", pt6, "ptE", ptE);

    if abs(frameRot) > 0
        R = [cos(frameRot), -sin(frameRot); ...
             sin(frameRot),  cos(frameRot)];
        pt0 = R*pt0; pt1 = R*pt1; pt2 = R*pt2; pt3 = R*pt3; pt4 = R*pt4;
        pt5 = R*pt5; pt6 = R*pt6; ptE = R*ptE;
    end

    vx = gradient(ptE(1,:), timeAxis);
    vy = gradient(ptE(2,:), timeAxis);
    speed = hypot(vx, vy);

    cycleState = struct();
    cycleState.timeAxis = timeAxis;
    cycleState.crankAngles = crankAngles;
    cycleState.gaitPercent = gaitPercent;
    cycleState.angleMat = angleMat;
    cycleState.qHist = qHist;
    cycleState.exitflag = exitflag;
    cycleState.resnorm = resnorm;
    cycleState.pt0 = pt0; cycleState.pt1 = pt1; cycleState.pt2 = pt2; cycleState.pt3 = pt3; cycleState.pt4 = pt4;
    cycleState.pt5 = pt5; cycleState.pt6 = pt6; cycleState.ptE = ptE;
    cycleState.rawMechanismFrame = raw;
    cycleState.frameRotation_rad = frameRot;
    cycleState.vx = vx; cycleState.vy = vy; cycleState.speed = speed;
    cycleState.span = [spanOf(ptE(1,:)); spanOf(ptE(2,:))];
    cycleState.area = abs(polyarea(ptE(1,:), ptE(2,:)));
end

function residualVec = loopEquations(unknownAngles, theta1, linkLen)
%LOOPEQUATIONS Five vector-loop closures, two scalar equations each.
    u = @(a) [cos(a); sin(a)];

    theta2 = unknownAngles(1); theta3 = unknownAngles(2); theta5 = unknownAngles(3); theta6 = unknownAngles(4);
    theta7 = unknownAngles(5); theta8 = unknownAngles(6); theta9 = unknownAngles(7); theta10 = unknownAngles(8);
    theta11 = unknownAngles(9); theta12 = unknownAngles(10);

    e1 = u(theta1);  e2 = u(theta2);  e3 = u(theta3);  e4 = [1;0];
    e5 = u(theta5);  e6 = u(theta6);  e7 = u(theta7);  e8 = u(theta8);
    e9 = u(theta9);  e10 = u(theta10); e11 = u(theta11); e12 = u(theta12);

    Fupper = linkLen(1)*e1 + linkLen(2)*e2 - linkLen(3)*e3 - linkLen(4)*e4;
    Flower = linkLen(1)*e1 + linkLen(7)*e7 - linkLen(8)*e8 - linkLen(4)*e4;
    Ftriangle = linkLen(3)*e3 + linkLen(5)*e5 - linkLen(6)*e6;
    Fpara = linkLen(8)*e8 + linkLen(9)*e9 - linkLen(6)*e6 - linkLen(10)*e10;
    Ffoot = linkLen(11)*e11 + linkLen(12)*e12 - linkLen(9)*e9;

    residualVec = [Fupper; Flower; Ftriangle; Fpara; Ffoot];
end

function jacobianMat = loopJacobian(unknownAngles, linkLen)
%LOOPJACOBIAN Analytic Jacobian dF/dq for the five vector-loop equations.
    d = @(a) [-sin(a); cos(a)];

    theta2 = unknownAngles(1); theta3 = unknownAngles(2); theta5 = unknownAngles(3); theta6 = unknownAngles(4);
    theta7 = unknownAngles(5); theta8 = unknownAngles(6); theta9 = unknownAngles(7); theta10 = unknownAngles(8);
    theta11 = unknownAngles(9); theta12 = unknownAngles(10);

    jacobianMat = zeros(10,10);

    % Upper loop: L1*e1 + L2*e2 - L3*e3 - L4*e4 = 0.
    jacobianMat(1:2,1) =  linkLen(2)*d(theta2);
    jacobianMat(1:2,2) = -linkLen(3)*d(theta3);

    % Lower loop: L1*e1 + L7*e7 - L8*e8 - L4*e4 = 0.
    jacobianMat(3:4,5) =  linkLen(7)*d(theta7);
    jacobianMat(3:4,6) = -linkLen(8)*d(theta8);

    % Coupler triangle: L3*e3 + L5*e5 - L6*e6 = 0.
    jacobianMat(5:6,2) =  linkLen(3)*d(theta3);
    jacobianMat(5:6,3) =  linkLen(5)*d(theta5);
    jacobianMat(5:6,4) = -linkLen(6)*d(theta6);

    % Parallelogram-like loop: L8*e8 + L9*e9 - L6*e6 - L10*e10 = 0.
    jacobianMat(7:8,4) = -linkLen(6)*d(theta6);
    jacobianMat(7:8,6) =  linkLen(8)*d(theta8);
    jacobianMat(7:8,7) =  linkLen(9)*d(theta9);
    jacobianMat(7:8,8) = -linkLen(10)*d(theta10);

    % Foot triangle: L11*e11 + L12*e12 - L9*e9 = 0.
    jacobianMat(9:10,7) = -linkLen(9)*d(theta9);
    jacobianMat(9:10,9) =  linkLen(11)*d(theta11);
    jacobianMat(9:10,10) = linkLen(12)*d(theta12);
end

function nonlinearSolver = makeLoopSolver()
%MAKELOOPSOLVER Use fsolve when available; otherwise use a built-in nonlinearSolver.
    nonlinearSolver = struct();
    nonlinearSolver.useFsolve = exist("fsolve", "file") == 2;
    nonlinearSolver.tolF = 1e-10;
    nonlinearSolver.tolStep = 1e-11;
    nonlinearSolver.maxIter = 80;
    nonlinearSolver.maxLineSearch = 16;

    if nonlinearSolver.useFsolve
        nonlinearSolver.options = optimoptions("fsolve", ...
            "Display", "off", ...
            "FunctionTolerance", nonlinearSolver.tolF, ...
            "StepTolerance", nonlinearSolver.tolStep, ...
            "OptimalityTolerance", nonlinearSolver.tolF, ...
            "MaxIterations", 120, ...
            "MaxFunctionEvaluations", 1500);
    else
        nonlinearSolver.options = [];
    end
end

function [x, residualVec, exitflag] = solveLoopAngles(fun, jac, x0, nonlinearSolver)
%SOLVELOOPANGLES Nonlinear solve wrapper.
    if nonlinearSolver.useFsolve
        [x, residualVec, exitflag] = fsolve(fun, x0, nonlinearSolver.options);
        return;
    end

    [x, residualVec, exitflag] = dampedNewtonSolve(fun, jac, x0, nonlinearSolver);
end

function [x, residualVec, exitflag] = dampedNewtonSolve(fun, jac, x0, nonlinearSolver)
%DAMPEDNEWTONSOLVE Small square-system Newton nonlinearSolver for this linkage.
% It is included so the script remains runnable without Optimization Toolbox.
    x = x0(:);
    residualVec = fun(x);
    nrm = norm(residualVec);
    exitflag = 0;

    for iter = 1:nonlinearSolver.maxIter
        if nrm < nonlinearSolver.tolF
            exitflag = 1;
            return;
        end

        jacobianMat = jac(x);
        step = -jacobianMat \ residualVec;
        if any(~isfinite(step)) || norm(step) > 5
            step = -pinv(jacobianMat) * residualVec;
        end

        alpha = 1.0;
        accepted = false;
        for ls = 1:nonlinearSolver.maxLineSearch
            xt = x + alpha*step;
            Ft = fun(xt);
            nt = norm(Ft);
            if nt < nrm || nt < nonlinearSolver.tolF
                x = xt;
                residualVec = Ft;
                nrm = nt;
                accepted = true;
                break;
            end
            alpha = 0.5*alpha;
        end

        if ~accepted
            x = x + alpha*step;
            residualVec = fun(x);
            nrm = norm(residualVec);
        end

        if norm(alpha*step) < nonlinearSolver.tolStep*(1 + norm(x))
            exitflag = double(nrm < 1e-7);
            return;
        end
    end
end

function pointSet = positionsFromAngles(linkLen, linkAngles)
%POSITIONSFROMANGLES Forward kinematics from solved link orientations.
    u = @(a) [cos(a); sin(a)];

    pt0 = [0;0];
    pt3 = [linkLen(4);0];
    pt1 = pt0 + linkLen(1)*u(linkAngles(1));
    pt2 = pt1 + linkLen(2)*u(linkAngles(2));
    pt5 = pt1 + linkLen(7)*u(linkAngles(7));
    pt4 = pt2 + linkLen(5)*u(linkAngles(5));
    pt6 = pt5 + linkLen(9)*u(linkAngles(9));
    ptE = pt5 + linkLen(11)*u(linkAngles(11));

    pointSet = struct("pt0", pt0, "pt1", pt1, "pt2", pt2, "pt3", pt3, ...
        "pt4", pt4, "pt5", pt5, "pt6", pt6, "ptE", ptE);
end

function unknownAngles = geometricInitialGuess(linkLen, theta1)
%GEOMETRICINITIALGUESS Build the physical assembly mode from circle intersections.
% This is not the nonlinearSolver; it selects the same branch before fsolve refines
% the vector-loop equations.
    pt0 = [0;0];
    pt3 = [linkLen(4);0];
    pt1 = pt0 + linkLen(1)*[cos(theta1); sin(theta1)];

    P2c = circleIntersections(pt1, linkLen(2), pt3, linkLen(3));
    pt2 = pickBy(P2c, "maxY");

    P5c = circleIntersections(pt1, linkLen(7), pt3, linkLen(8));
    pt5 = pickBy(P5c, "minY");

    P4c = circleIntersections(pt2, linkLen(5), pt3, linkLen(6));
    pt4 = pickBy(P4c, "maxX");

    P6c = circleIntersections(pt5, linkLen(9), pt4, linkLen(10));
    pt6 = pickBy(P6c, "maxX");
    if pt6(2) > min(pt4(2), pt5(2))
        pt6 = pickBy(P6c, "minY");
    end

    PEc = circleIntersections(pt5, linkLen(11), pt6, linkLen(12));
    ptE = pickBy(PEc, "minY");

    ang = @(a,b) atan2(b(2)-a(2), b(1)-a(1));
    unknownAngles = [
        ang(pt1, pt2)
        ang(pt3, pt2)
        ang(pt2, pt4)
        ang(pt3, pt4)
        ang(pt1, pt5)
        ang(pt3, pt5)
        ang(pt5, pt6)
        ang(pt4, pt6)
        ang(pt5, ptE)
        ang(ptE, pt6)
    ];
end

function pts = circleIntersections(c1, r1, c2, r2)
%CIRCLEINTERSECTIONS Return the two intersection points of two circles.
    dvec = c2 - c1;
    d = norm(dvec);
    if d < eps
        error("Coincident circle centers in initial guess.");
    end
    if d > r1 + r2 || d < abs(r1 - r2)
        error("The selected link lengths cannot assemble at this crankAngles angle.");
    end

    a = (r1^2 - r2^2 + d^2) / (2*d);
    h2 = max(r1^2 - a^2, 0);
    h = sqrt(h2);
    ex = dvec / d;
    ey = [-ex(2); ex(1)];
    p = c1 + a*ex;
    pts = [p + h*ey, p - h*ey];
end

function p = pickBy(pts, mode)
%PICKBY Select one of two branch points by geometric criterion.
    switch mode
        case "maxY"
            [~, idx] = max(pts(2,:));
        case "minY"
            [~, idx] = min(pts(2,:));
        case "maxX"
            [~, idx] = max(pts(1,:));
        case "minX"
            [~, idx] = min(pts(1,:));
        otherwise
            error("Unknown branch-selection mode.");
    end
    p = pts(:,idx);
end

function unknownAngles = unwrapNear(qNew, qOld)
%UNWRAPNEAR Keep periodic angle solutions close to the previous step.
    if isempty(qOld) || any(~isfinite(qOld))
        unknownAngles = qNew;
        return;
    end
    unknownAngles = qNew + 2*pi*round((qOld - qNew)/(2*pi));
end

function refCurve = makeReferenceGait(gaitPercent, span)
%MAKEREFERENCEGAIT Smooth gait-like reference with the paper-reported spans.
% This is for visual validationGrid only; it is not the unpublished human dataset.
    s = gaitPercent(:).' / 100;
    knot = [0.00 0.08 0.16 0.28 0.42 0.58 0.72 0.86 1.00];
    xShape = [0.48 0.42 0.25 -0.15 -0.50 -0.35 -0.10 0.20 0.48];
    yShape = [0.55 0.95 1.00 0.35 0.15 0.06 0.00 0.18 0.55];
    x = pchip(knot, xShape, s);
    y = pchip(knot, yShape, s);

    x = span(1) * (x - min(x)) / spanOf(x);
    y = span(2) * (y - min(y)) / spanOf(y);
    x = x - mean(x);
    y = y - min(y);

    refCurve = [x; y];
end

function refAligned = alignCurveToSimulation(refCurve, ptE)
%ALIGNCURVETOSIMULATION Translate reference to same lower-left envelope.
    refAligned = refCurve;
    refAligned(1,:) = refCurve(1,:) - mean(refCurve(1,:)) + mean(ptE(1,:));
    refAligned(2,:) = refCurve(2,:) - min(refCurve(2,:)) + min(ptE(2,:));
end

function plotPaperStyleFigure4(trajData, refCurve, linkLen)
%PLOTPAPERSTYLEFIGURE4 Match the composite mechanism/gait plot in the paper.
    figH = figure("Name", "Paper-style Fig. 4: mechanism and gait curves", ...
        "Color", "w", "Position", [80 80 1120 560]);

    axMech = axes(figH, "Position", [0.06 0.16 0.53 0.76]);
    axes(axMech);
    k = 22;
    drawPaperMechanism(trajData, k, linkLen);
    hold on;

    sampleIdx = 1:5:numel(trajData.gaitPercent);
    plot(refCurve(1,sampleIdx), refCurve(2,sampleIdx), "cx", ...
        "MarkerSize", 5.5, "LineWidth", 1.1);
    plot(trajData.ptE(1,sampleIdx), trajData.ptE(2,sampleIdx), "ko", ...
        "MarkerSize", 4.0, "LineWidth", 1.1);
    text(trajData.ptE(1,end)+1.0, trajData.ptE(2,end), "Endpoint", ...
        "FontSize", 10, "VerticalAlignment", "middle");

    linkText = makeLengthListText(linkLen.val);
    txtX = min([trajData.pt0(1,:), trajData.ptE(1,:)]) - 6;
    txtY = max([trajData.pt2(2,:), trajData.pt4(2,:)]) - 4;
    text(txtX, txtY, linkText, "FontSize", 9, "FontName", "Consolas", ...
        "VerticalAlignment", "top", "Interpreter", "none");

    axis equal;
    axis off;
    title("Parameterized 12-link modified Jansen mechanism", ...
        "FontSize", 11, "FontWeight", "normal");

    caption = sprintf(['Fig. 4  Optimized to the gait envelope, the structure produces an endpoint trajectory with ' ...
        'x-span %.2f cm and y-span %.2f cm. Crosses denote the reference envelope; circles/solid curves denote the ' ...
        'simulated mechanism trajectory under constant crankAngles speed.'], trajData.span(1), trajData.span(2));
    annotation(figH, "textbox", [0.08 0.015 0.86 0.10], ...
        "String", caption, ...
        "EdgeColor", "none", "FontSize", 10, "FontWeight", "bold");

    % Separate window for gait-cycle plots
    fig2 = figure("Name", "Gait Cycle Curves", ...
        "Color", "w", "Position", [850 100 520 620]); %#ok<NASGU>

    tiledlayout(2,1,"Padding","compact","TileSpacing","compact");

    % X-axis graph
    nexttile;

    predX = trajData.ptE(1,:) - mean(trajData.ptE(1,:));
    refX = refCurve(1,:) - mean(refCurve(1,:));

    plot(trajData.gaitPercent, refX, "g--", "LineWidth", 1.5); hold on;
    plot(trajData.gaitPercent, predX, "c-", "LineWidth", 1.8);
    grid on;
    xlim([0 100]);
    ylim([-40 40]);
    ylabel("x-axis (cm)");
    title("Horizontal Motion Over Gait Cycle");
    set(gca, "FontSize", 9);

    % Y-axis graph
    nexttile;

    predY = trajData.ptE(2,:) - min(trajData.ptE(2,:));
    refY = refCurve(2,:) - min(refCurve(2,:));

    plot(trajData.gaitPercent, refY, "g--", "LineWidth", 1.5); hold on;
    plot(trajData.gaitPercent, predY, "c-", "LineWidth", 1.8);
    grid on;
    xlim([0 100]);
    ylim([0 15]);
    xlabel("Gait Cycle (%)");
    ylabel("y-axis (cm)");
    title("Vertical Motion Over Gait Cycle");
    legend("meta-trajectory", "predicted trajectory", ...
        "Location", "northeast", "FontSize", 8);
    set(gca, "FontSize", 9);

    axes(axMech);
end

function drawPaperMechanism(trajData, k, linkLen)
%DRAWPAPERMECHANISM Draw one pose with labels similar to the paper figure.
    linePairs = {
        "pt0","pt1","L_1"; "pt1","pt2","L_2"; "pt3","pt2","L_3"; "pt0","pt3","L_4";
        "pt2","pt4","L_5"; "pt3","pt4","L_6"; "pt1","pt5","L_7"; "pt3","pt5","L_8";
        "pt5","pt6","L_9"; "pt4","pt6","L_{10}"; "pt5","ptE","L_{11}"; "ptE","pt6","L_{12}"};

    for i = 1:size(linePairs,1)
        A = trajData.(linePairs{i,1})(:,k);
        B = trajData.(linePairs{i,2})(:,k);
        plot([A(1), B(1)], [A(2), B(2)], "k-", "LineWidth", 1.5);
        hold on;
        mid = 0.52*A + 0.48*B;
        text(mid(1), mid(2), linePairs{i,3}, ...
            "FontSize", 11, "FontWeight", "bold", "Interpreter", "tex");
    end

    pointNames = ["pt0","pt1","pt2","pt3","pt4","pt5","pt6","ptE"];
    P = zeros(2, numel(pointNames));
    for i = 1:numel(pointNames)
        P(:,i) = trajData.(pointNames(i))(:,k);
    end
    plot(P(1,:), P(2,:), "ko", "MarkerFaceColor", "c", "MarkerSize", 4.5);

    xMargin = 8;
    yMargin = 8;
    xlim([min([P(1,:), trajData.ptE(1,:)])-xMargin, max([P(1,:), trajData.ptE(1,:)])+xMargin]);
    ylim([min([P(2,:), trajData.ptE(2,:)])-yMargin, max([P(2,:), trajData.ptE(2,:)])+yMargin]);

    % Put link labels a little away from the length-list text.
    unused = linkLen; %#ok<NASGU>
end

function txt = makeLengthListText(linkLen)
%MAKELENGTHLISTTEXT Link length annotation used in the paper-style figure.
    lines = strings(14,1);
    lines(1) = "Unit: cm";
    for i = 1:12
        lines(i+1) = sprintf("linkLen%-2d = %4.1f", i, linkLen(i));
    end
    txt = strjoin(lines, newline);
end

function validationGrid = plotNinePatternValidation(Lnom, cfg)
%PLOTNINEPATTERNVALIDATION Create a 3x3 paper-style RMSE validationGrid grid.
% The PDFs do not contain the original 113-subject ankle database. This grid
% uses nine gait-envelope variants and solves the actual mechanism for each.
    xSpanGrid = [56 53 50; 53 50 47; 50 47 44];
    ySpanGrid = [14.2 13.6 13.0; 13.6 13.0 12.4; 13.0 12.4 11.8];
    L1Grid = Lnom(1) + [0.45 0.25 0.05; 0.25 0.00 -0.20; 0.05 -0.20 -0.45];
    L4Grid = Lnom(4) + [-0.90 -0.50 -0.10; -0.50 0.00 0.50; -0.10 0.50 0.90];
    L8Grid = Lnom(8) + [0.90 0.50 0.10; 0.50 0.00 -0.50; 0.10 -0.50 -0.90];

    figH = figure("Name", "Paper-style 3x3 validationGrid grid", ...
        "Color", "w", "Position", [120 90 1040 620]);
    tileLayout = tiledlayout(3,3, "Padding", "compact", "TileSpacing", "compact");

    validationGrid = struct();
    validationGrid.patterns = repmat(struct("linkLen", [], "reference", [], ...
        "simulation", [], "rmse", [], "span", []), 3, 3);

    for row = 1:3
        for col = 1:3
            linkLen = Lnom;
            linkLen(1) = L1Grid(row,col);
            linkLen(4) = L4Grid(row,col);
            linkLen(8) = L8Grid(row,col);
            cycleState = simulateJansenCycle(linkLen, cfg.NAdvanced, cfg.omega, cfg.phaseOffset, [], ...
                cfg.gaitFrameRotationDeg*pi/180);

            refCurve = makeReferenceGait(cycleState.gaitPercent, [xSpanGrid(row,col); ySpanGrid(row,col)]);
            refCurve = alignCurveToSimulation(refCurve, cycleState.ptE);

            [simPanelCurve, refPanelCurve] = panelNormalizeCurves(cycleState.ptE, refCurve);
            rmse = sqrt(mean(sum((simPanelCurve - refPanelCurve).^2, 1)));

            axH = nexttile;
            sampleIdx = 1:2:numel(cycleState.gaitPercent);
            plot(refPanelCurve(1,sampleIdx), refPanelCurve(2,sampleIdx), "cx", ...
                "MarkerSize", 5.0, "LineWidth", 1.0); hold on;
            plot(simPanelCurve(1,sampleIdx), simPanelCurve(2,sampleIdx), "ko", ...
                "MarkerSize", 4.0, "LineWidth", 1.0);
            grid on;
            xlim([0 70]);
            ylim([0 20]);
            text(5, 17, sprintf("RMSE=%.2f", rmse), ...
                "FontSize", 11, "FontWeight", "normal");
            set(axH, "FontSize", 8);

            if row < 3
                set(axH, "XTickLabel", []);
            else
                xlabel("x- axis (cm)");
            end
            if col > 1
                set(axH, "YTickLabel", []);
            else
                ylabel("y- axis (cm)");
            end
            if row == 1 && col == 1
                legend("reference", "simulation", "Location", "northwest", "FontSize", 7);
            end

            validationGrid.patterns(row,col).linkLen = linkLen;
            validationGrid.patterns(row,col).reference = refPanelCurve;
            validationGrid.patterns(row,col).simulation = simPanelCurve;
            validationGrid.patterns(row,col).rmse = rmse;
            validationGrid.patterns(row,col).span = cycleState.span;
        end
    end

    title(tileLayout, "Reference-vs-simulation endpoint trajectories for nine gait envelopes", ...
        "FontSize", 12, "FontWeight", "bold");
end

function [simPanelCurve, refPanelCurve] = panelNormalizeCurves(simCurve, refCurve)
%PANELNORMALIZECURVES Translate curves to a common 0-to-70 cm plotting frame.
    simPanelCurve = simCurve;
    refPanelCurve = refCurve;

    xmin = min([simPanelCurve(1,:), refPanelCurve(1,:)]);
    ymin = min([simPanelCurve(2,:), refPanelCurve(2,:)]);

    simPanelCurve(1,:) = simPanelCurve(1,:) - xmin + 8;
    refPanelCurve(1,:) = refPanelCurve(1,:) - xmin + 8;
    simPanelCurve(2,:) = simPanelCurve(2,:) - ymin + 1.5;
    refPanelCurve(2,:) = refPanelCurve(2,:) - ymin + 1.5;
end

function animateJansenMechanism(trajData, linkLen, frameStep)
%ANIMATEJANSENMECHANISM Visible one-cycle mechanism simulation.
    if nargin < 3
        frameStep = 3;
    end

    if ~usejava("desktop")
        disp("Live animation skipped because MATLAB is running without the desktop UI.");
        return;
    end

    figH = figure("Name", "Running simulation: modified Jansen mechanism", ...
        "Color", "w", "Position", [160 120 900 560]);
    axH = axes("Parent", figH);

    allX = [trajData.pt0(1,:), trajData.pt1(1,:), trajData.pt2(1,:), trajData.pt3(1,:), trajData.pt4(1,:), ...
        trajData.pt5(1,:), trajData.pt6(1,:), trajData.ptE(1,:)];
    allY = [trajData.pt0(2,:), trajData.pt1(2,:), trajData.pt2(2,:), trajData.pt3(2,:), trajData.pt4(2,:), ...
        trajData.pt5(2,:), trajData.pt6(2,:), trajData.ptE(2,:)];
    xLim = [min(allX)-8, max(allX)+8];
    yLim = [min(allY)-8, max(allY)+8];

    for k = 1:frameStep:numel(trajData.gaitPercent)
        if ~isvalid(figH)
            break;
        end
        if ~isvalid(axH)
            axH = axes("Parent", figH);
        end
        cla(axH);
        axes(axH); %#ok<LAXES>
        drawMechanism(trajData, k, [0.80 0.35 0.10]);
        hold on;
        plot(trajData.ptE(1,1:k), trajData.ptE(2,1:k), "c-", "LineWidth", 2.2);
        plot(trajData.ptE(1,k), trajData.ptE(2,k), "ko", "MarkerFaceColor", "c", "MarkerSize", 7);
        grid on; axis equal;
        xlim(xLim); ylim(yLim);
        xlabel("x position (cm)");
        ylabel("y position (cm)");
        title(sprintf("Running one-DOF simulation, crankAngles angle = %.1f deg", ...
            trajData.crankAngles(k)*180/pi));
        text(xLim(1)+2, yLim(2)-4, ...
            sprintf("L1=%.2f cm, L4=%.2f cm, L8=%.2f cm", linkLen.val(1), linkLen.val(4), linkLen.val(8)), ...
            "BackgroundColor", "w", "Margin", 4);
        drawnow;
        pause(0.005);
    end
end

function plotEndEffector(trajData, refCurve, targetSpan, linkLen)
    figure("Name", "End-effector gait trajectory", "Color", "w");
    plot(trajData.ptE(1,:), trajData.ptE(2,:), "m-", "LineWidth", 2.4); hold on;
    plot(refCurve(1,:), refCurve(2,:), "g--", "LineWidth", 1.7);
    plot(trajData.ptE(1,1), trajData.ptE(2,1), "ko", "MarkerFaceColor", "k", "MarkerSize", 5);
    grid on; axis equal;
    xlabel("x position (cm)");
    ylabel("y position (cm)");
    title("Modified Jansen end-effector / ankle trajectory");
    legend("simulated mechanism", ...
        sprintf("gait-like reference %.2f x %.2f cm", targetSpan(1), targetSpan(2)), ...
        "start", "Location", "best");

    text(min(trajData.ptE(1,:)), max(trajData.ptE(2,:)), ...
        sprintf("L1=%.2f, L4=%.2f, L8=%.2f cm", linkLen.val(1), linkLen.val(4), linkLen.val(8)), ...
        "VerticalAlignment", "top", "BackgroundColor", "w", "Margin", 4);
end

function plotGaitCurves(trajData, refCurve)
    figure("Name", "Gait-cycle position curves", "Color", "w");
    tiledlayout(2,1, "Padding", "compact", "TileSpacing", "compact");

    nexttile;
    plot(trajData.gaitPercent, trajData.ptE(1,:), "m-", "LineWidth", 2.2); hold on;
    plot(trajData.gaitPercent, refCurve(1,:), "g--", "LineWidth", 1.5);
    grid on;
    ylabel("x (cm)");
    title("Horizontal ankle motion over gait cycle");
    legend("simulated", "reference envelope", "Location", "best");

    nexttile;
    plot(trajData.gaitPercent, trajData.ptE(2,:), "m-", "LineWidth", 2.2); hold on;
    plot(trajData.gaitPercent, refCurve(2,:), "g--", "LineWidth", 1.5);
    grid on;
    xlabel("gait cycle (%)");
    ylabel("y (cm)");
    title("Vertical ankle motion over gait cycle");
end

function plotVelocityCurves(trajData)
    figure("Name", "End-effector velocity curves", "Color", "w");
    tiledlayout(3,1, "Padding", "compact", "TileSpacing", "compact");

    nexttile;
    plot(trajData.gaitPercent, trajData.vx, "LineWidth", 1.8);
    grid on; ylabel("vx (cm/s)");
    title("End-effector velocity with constant crankAngles speed");

    nexttile;
    plot(trajData.gaitPercent, trajData.vy, "LineWidth", 1.8);
    grid on; ylabel("vy (cm/s)");

    nexttile;
    plot(trajData.gaitPercent, trajData.speed, "LineWidth", 1.8);
    grid on; xlabel("gait cycle (%)"); ylabel("|v| (cm/s)");
end

function plotConstraintResiduals(trajData)
    figure("Name", "Closed-loop validationGrid", "Color", "w");
    semilogy(trajData.gaitPercent, trajData.resnorm + eps, "k-", "LineWidth", 1.8);
    grid on;
    xlabel("gait cycle (%)");
    ylabel("||loop residual||_2 (cm)");
    title("Vector-loop closure residual from fsolve");
end

function plotMechanismSnapshots(trajData, linkLen)
    figure("Name", "Mechanism snapshots", "Color", "w");
    idx = round(linspace(1, numel(trajData.gaitPercent)-1, 7));
    for k = idx
        drawMechanism(trajData, k, [0.45 0.45 0.45]);
        hold on;
    end
    drawMechanism(trajData, idx(2), [0.20 0.55 0.20]);
    plot(trajData.ptE(1,:), trajData.ptE(2,:), "c-", "LineWidth", 2.0);
    grid on; axis equal;
    xlabel("x (cm)"); ylabel("y (cm)");
    title(sprintf("Assembly snapshots, linkLen = [%.2f %.1f %.1f %.2f ...] cm", ...
        linkLen.val(1), linkLen.val(2), linkLen.val(3), linkLen.val(4)));
end

function drawMechanism(trajData, k, color)
    linePairs = {
        "pt0","pt1"; "pt1","pt2"; "pt3","pt2"; "pt0","pt3";
        "pt1","pt5"; "pt3","pt5"; "pt2","pt4"; "pt3","pt4";
        "pt5","pt6"; "pt4","pt6"; "pt5","ptE"; "ptE","pt6"};

    for i = 1:size(linePairs,1)
        A = trajData.(linePairs{i,1})(:,k);
        B = trajData.(linePairs{i,2})(:,k);
        plot([A(1), B(1)], [A(2), B(2)], "-", "Color", color, "LineWidth", 1.2);
        hold on;
    end
    pts = [trajData.pt0(:,k), trajData.pt1(:,k), trajData.pt2(:,k), trajData.pt3(:,k), ...
        trajData.pt4(:,k), trajData.pt5(:,k), trajData.pt6(:,k), trajData.ptE(:,k)];
    plot(pts(1,:), pts(2,:), "o", "Color", color, "MarkerFaceColor", "w", "MarkerSize", 4);
end

function sensitivityData = runAdjustableLinkSensitivity(Lnom, cfg)
%RUNADJUSTABLELINKSENSITIVITY Quantify effects of L1, L4, and L8.
    adjustable = [1 4 8];
    pct = linspace(-0.04, 0.04, 9);
    spans = nan(numel(adjustable), numel(pct), 2);
    areas = nan(numel(adjustable), numel(pct));

    for i = 1:numel(adjustable)
        for j = 1:numel(pct)
            linkLen = Lnom;
            linkLen(adjustable(i)) = Lnom(adjustable(i)) * (1 + pct(j));
            cycleState = simulateJansenCycle(linkLen, cfg.NAdvanced, cfg.omega, cfg.phaseOffset, [], ...
                cfg.gaitFrameRotationDeg*pi/180);
            spans(i,j,:) = cycleState.span;
            areas(i,j) = cycleState.area;
        end
    end

    sensitivityData = struct();
    sensitivityData.adjustable = adjustable;
    sensitivityData.percentChange = 100*pct;
    sensitivityData.spans = spans;
    sensitivityData.areas = areas;
end

function plotSensitivity(sensitivityData)
    figure("Name", "Adjustable-link sensitivityData", "Color", "w");
    tiledlayout(1,2, "Padding", "compact", "TileSpacing", "compact");
    names = ["L1", "L4", "L8"];

    nexttile;
    for i = 1:3
        plot(sensitivityData.percentChange, squeeze(sensitivityData.spans(i,:,1)), ...
            "o-", "LineWidth", 1.6); hold on;
    end
    grid on; xlabel("link length change (%)"); ylabel("x-span / stride (cm)");
    title("Effect on horizontal span");
    legend(names, "Location", "best");

    nexttile;
    for i = 1:3
        plot(sensitivityData.percentChange, squeeze(sensitivityData.spans(i,:,2)), ...
            "o-", "LineWidth", 1.6); hold on;
    end
    grid on; xlabel("link length change (%)"); ylabel("y-span / step height (cm)");
    title("Effect on vertical span");
    legend(names, "Location", "best");
end

function lsMap = runLeastSquaresSpanMap(Lnom, cfg, targetSpan)
%RUNLEASTSQUARESSPANMAP Demonstrates Lambda = Psi*Sigma from the paper.
% Lambda stores [L1; L4; L8] samples. Sigma stores [xspan; yspan] samples.
    L1set = Lnom(1) + [-0.6 0 0.6];
    L4set = Lnom(4) + [-1.0 0 1.0];
    L8set = Lnom(8) + [-1.0 0 1.0];

    n = numel(L1set) * numel(L4set) * numel(L8set);
    Lambda = nan(3, n);
    Sigma = nan(2, n);
    samples = nan(n, 5);
    c = 0;

    for a = 1:numel(L1set)
        for b = 1:numel(L4set)
            for d = 1:numel(L8set)
                c = c + 1;
                linkLen = Lnom;
                linkLen(1) = L1set(a);
                linkLen(4) = L4set(b);
                linkLen(8) = L8set(d);
                cycleState = simulateJansenCycle(linkLen, cfg.NAdvanced, cfg.omega, cfg.phaseOffset, [], ...
                    cfg.gaitFrameRotationDeg*pi/180);
                Lambda(:,c) = [linkLen(1); linkLen(4); linkLen(8)];
                Sigma(:,c) = cycleState.span;
                samples(c,:) = [linkLen(1), linkLen(4), linkLen(8), cycleState.span(1), cycleState.span(2)];
            end
        end
    end

    Psi = Lambda * pinv(Sigma);
    LpredAdj = Psi * targetSpan;
    Lpred = Lnom;
    Lpred([1 4 8]) = LpredAdj;
    predicted = simulateJansenCycle(Lpred, cfg.N, cfg.omega, cfg.phaseOffset, [], ...
        cfg.gaitFrameRotationDeg*pi/180);

    lsMap = struct();
    lsMap.Lambda = Lambda;
    lsMap.Sigma = Sigma;
    lsMap.Psi = Psi;
    lsMap.samples = samples;
    lsMap.Lpred = Lpred;
    lsMap.predicted = predicted;
end

function plotLeastSquaresDemo(lsMap, targetSpan)
    figure("Name", "Least-squares span-to-link map", "Color", "w");
    tiledlayout(1,2, "Padding", "compact", "TileSpacing", "compact");

    nexttile;
    scatter(lsMap.Sigma(1,:), lsMap.Sigma(2,:), 45, "k", "filled"); hold on;
    plot(targetSpan(1), targetSpan(2), "mp", "MarkerSize", 14, "MarkerFaceColor", "m");
    plot(lsMap.predicted.span(1), lsMap.predicted.span(2), "co", "MarkerSize", 9, "LineWidth", 2);
    grid on;
    xlabel("x-span (cm)");
    ylabel("y-span (cm)");
    title("Sampled span cloud and requested span");
    legend("simulation samples", "requested span", "span from predicted links", "Location", "best");

    nexttile;
    plot(lsMap.predicted.ptE(1,:), lsMap.predicted.ptE(2,:), "m-", "LineWidth", 2.2);
    grid on; axis equal;
    xlabel("x (cm)"); ylabel("y (cm)");
    title(sprintf("LS predicted L1=%.2f, L4=%.2f, L8=%.2f cm", ...
        lsMap.Lpred(1), lsMap.Lpred(4), lsMap.Lpred(8)));
end

function s = spanOf(v)
%SPANOF Max-min span without relying on toolbox-specific range behavior.
    s = max(v(:)) - min(v(:));
end
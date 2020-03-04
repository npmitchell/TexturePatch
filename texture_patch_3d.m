function texture_patch_3d( FF, VV, TF, TV, IV, Options)
%TEXTURE_PATCH_3D This function will display a triangulated mesh, like the
%MATLAB function "patch", with a texture mapping given by the 3D image
%volume I. One set of inputs to
%this function is the "physical" or "real-space" mesh triangulation. This
%is the mesh triangulation you actually wish to display.  The other input
%is the "virtual" or "texture-space" mesh triangulation.  This input maps
%the physical mesh into texture space.  The virtual mesh must have the same
%number of facets as the physical mesh (each face needs a texture
%mapping!), but may have different connectivities (i.e. the physical mesh
%may be simply connected while the virtual mesh is multiply connected) or a
%different number of vertices (i.e. re-using the texturing from a single
%template face for multiple physical faces).
%
%   texture_patch_3d( FF, VV, TF, TV, IV, Options );
%
%   Input Parameters:
%       - FF:       #Fx3 face connectivity list of the physical mesh
%       - VV:       #Vx3 vertex coordinate list of the physical mesh. NOTE:
%                   VV is provided in (X,Y,Z) format!
%       - TF:       #Fx3 face connectivity list of the virtual mesh
%       - TV:       #Kx3  texture coordinate list. Coordinates can be given
%                   as real pixel positions in the texture image or as a
%                   range [0..1]. NOTE: TV is provided in
%                   (row, column, page) format!
%       - IV:        The 3D texture image volume.This function supports:
%                   Grayscale:  [I1 x I2 x I3]
%                   RGB:    { [I1 x I2 x I3] } x 3 cell array (3D)
%
%       - Options:  Structure containing the standard options for a
%                   textured surface patch, such as EdgeColor, EdgeAlpha,
%                   etc.  See MATLAB documentation for more information.
%
%       - Options.PSize:    Special option, defines the image texture size
%       for each individual polygon.  A lower number gives a more 
%       pixelated texture {64}
%
%       - Options.Rotation: Special option, defines the rotation matrix of 
%       all surface patches
%
%       - Options.Translation: Special option, defines the translation 
%       vector applied after rotation on all surfaces
%
%       - Options.Rotation: Special option, defines the dilation factor 
%       applied to all surfaces after rotation and translation
%
%       - Options.ApplyAmbientOcclusion: Determines if the texture
%       colors should be modified by the ambient occlusion of the 
%       underlying triangulation {'false'}
%
%       - Options.AmbientOcclusion:    #Vx1 list of ambient occlusion
%       values
%
%       - Options.AmbientOcculsionFactor: A scalar between [0,1]
%       determining how strongly the ambient occlusion modifies the 
%       texture mapping {1}
%
%       - Options.AmbientOcculsionSamples: The number of samples to use
%       when computing the ambient occlusion {1000}
%
%       - Options.Unoriented: Treat the surface as unoriented when
%       applying ambient occlusion {'false'}
%
%       - Options.AddLights: Will add lighting to the plot {'false'}
%
%       - Options.ScalarField: A vertex-based or face based scalar field
%       used to augment texture mapping with external color
%
%       - Options.ScalarCLim: The colormap limits used to display the
%       scalar field. By default it will use the full range of values in
%       the scalar field
%
%       - Options.ScalarColorMap: The colormap used to display the scalar
%       field.  Either a a string corresponding to a built-in MATLAB
%       colormap or a user supplied #Cx3 colormap {'parula'}
%
%       - Options.ScalarAlpha: A scalar between [0,1] determining
%       how the transparency of the overlay of the scalar field mapping
%       {0.1};
%
%       - Options.VertexNormals: #Vx3 list of vertex unit normals
%
% See also
% --------
% texture_patch_to_image.m
%
%   by Dillon Cislo 08/14/2019 & NPMichell 8/19-2/20
%   NPMitchell added Rotation, Translation, & Dilation options 09/2019
%   NPMitchell grouped surfaces into a Parent container for speedup 09/2019
%   NPMitchell added colorize to options 12/2019
%   NPMitchell added capability for RGB color input images
%   NPMitchell added Imax, Imin for clipping the interpolated intensities

%--------------------------------------------------------------------------
% INPUT PROCESSING
%--------------------------------------------------------------------------

% Validate input triangulations -------------------------------------------

validateattributes(FF, {'numeric'}, ...
    {'2d', 'ncols', 3, 'finite', 'nonnan', 'real', 'integer', 'positive'});
validateattributes(VV, {'numeric'}, ...
    {'2d', 'finite', 'nonnan', 'real'});
validateattributes(TF, {'numeric'}, ...
    {'2d', 'ncols', 3, 'finite', 'nonnan', 'real', 'integer', 'positive'});
validateattributes(TV, {'numeric'}, ...
    {'2d', 'finite', 'nonnan', 'real'});

% Check that the number of faces is consistent
if ( size(FF,2) ~= size(TF,2) )
    error('texture_patch_3d:inputs', ...
        [ 'Number of real-space faces must match', ...
        'the number of texture-space faces' ] );
end

% Check the dimensions of the physical vertices
if ( size(VV,2) ~= 3 )
    if( size(VV,2) == 2 )
        VV = [ VV zeros( size(VV,1), 1 ) ];
    else
        error('texture_patch_3d:inputs', ...
            'Invalid real-space vertex input');
    end
end

% Check the dimensions of the texture vertices
if ( size(TV,2) ~= 3 )
    error('texture_patch_3d:inputs', 'Invalid texture vertex input');
end

% Re-scale texture vertices to true pixel positions if necessary
if ( max(TV(:)) < 2 )
    
    TV = [ (size(IV,1)-1) .* TV(:,1) + 1, ...
        (size(IV,2)-1) .* TV(:,2) + 1, ...
        (size(IV,3)-1) .* TV(:,2) + 1 ];
    
end

% Validate texture patch options ------------------------------------------
if (nargin < 6), Options = Struct(); end

% Patch FaceColor MUST be a texture.
Options.FaceColor = 'texturemap';

% Size of the 2D texture image used for each triangle
if isfield( Options, 'PSize' )
    sizep = round(Options.PSize(1));
    Options = rmfield(Options, 'PSize');
else
    sizep = 64;
end

% Global rotation for all faces
if isfield( Options, 'Rotation' )
    rot = Options.Rotation ;
    Options = rmfield(Options, 'Rotation');
    rotate = true ;
else
    rotate = false ;
end

% Global translation for all faces
if isfield( Options, 'Translation' )
    trans = Options.Translation ;
    Options = rmfield(Options, 'Translation');
    translate = true ;
else
    translate = false ;
end

% Global dilation for all faces
if isfield( Options, 'Dilation' )
    dilation = Options.Dilation ;
    Options = rmfield( Options, 'Dilation' );
    dilate = true ;
else
    dilate = false ;
end

% Determine if input is RGB or grayscale
if isfield( Options, 'isRGB' )
    isRGB = Options.isRGB;
    Options = rmfield( Options, 'isRGB' );
else
    if ( iscell(IV) || ( (size(TV,2) == 2) && (size(IV,3) == 3) ) )
        isRGB = true;
    else
        isRGB = false;
    end
end

% Determine if input is RGB or grayscale
if isfield( Options, 'Imax' )
    Imax = Options.Imax;
    Options = rmfield( Options, 'Imax' );
else
    Imax = Inf ;
end
if isfield( Options, 'Imin' )
    Imin = Options.Imin;
    Options = rmfield( Options, 'Imin' );
else
    Imin = -Inf ;
end

% Apply scaling of grayscape texture data by scalar field
% if isRGB
%     colorize = false ;
% else
%     if isfield( Options, 'FaceScalarField' )
%         heat = Options.FaceScalarField ;
%         Options = rmfield( Options, 'FaceScalarField' );
%         colorize = true ;
%     else
%         colorize = false ;
%     end
% end

% A vertex- or face-based scalar field used to colorize data
if isfield( Options, 'ScalarField' )
    S = Options.ScalarField;
    Options = rmfield( Options, 'ScalarField' );
    colorize = true;
    
    if ~( (numel(S) == size(VV,1)) || ...
            (numel(S) == size(FF,1)) )
        error('texture_patch_3d:inputs', ...
            'Invalid scalar field colorization input');
    end
    
else
    colorize = false;
    S = [];
end

% Limits of the scalar field colormap
if isfield( Options, 'ScalarCLim' )
    SCLim = Options.ScalarCLim;
    Options = rmfield( Options, 'ScalarCLim' );
else
    if colorize
        SCLim = [ min(S) max(S) ];
    end
end

% The scalar field colormap
if isfield( Options, 'ScalarColorMap' )
    SCMap = Options.ScalarColorMap;
    Options = rmfield( Options, 'ScalarColorMap' );
    
    if isnumeric(SCMap)
        if ~(ismatrix(SCMap) && (size(SCMap,2) == 3))
            error('texture_patch_3d:inputs', ...
                'Invalid scalar field colormap input');
        end
    else
        if ~exist('SCMap', 'var')
            error('texture_patch_3d:inputs', ...
                'Invalid scalar field colormap input');
        end
        SCMap = eval([SCMap '(256)']);
    end
    
else
    SCMap = parula; % 256 x 3
end

% Weights for combining the scalar field with the texture mapping
if isfield( Options, 'ScalarAlpha' )
    SAlpha = Options.ScalarAlpha;
    Options = rmfield( Options, 'ScalarAlpha' );
else
    SAlpha = 0.1;
end

% Apply ambient occlusion
if isfield( Options, 'ApplyAmbientOcclusion' )
    applyAO = Options.ApplyAmbientOcclusion;
    Options =rmfield(Options, 'ApplyAmbientOcclusion');
else
    applyAO = false;
end

% Ambient occlusion values
if isfield( Options, 'AmbientOcclusion' )
    AO = Options.AmbientOcclusion;
    Options = rmfield( Options, 'AmbientOcclusion' );
    if ~isempty(AO), applyAO = true; end
else
    AO = [];
end

% Ambient occlusion factor
if isfield( Options, 'AmbientOcclusionFactor' )
    AOFactor = Options.AmbientOcclusionFactor;
    Options = rmfield( Options, 'AmbientOcclusionFactor' );
else
    AOFactor = 1;
end

% Ambient occlusion sample number
if isfield( Options, 'AmbientOcclusionSamples' )
    AOSamples = Options.AmbientOcclusionSamples;
    Options = rmfield( Options, 'AmbientOcclusionSamples' );
else
    AOSamples = 1000;
end

% Surface orientation handling
if isfield( Options, 'Unoriented' )
    unoriented = Options.Unoriented;
    Options = rmfield( Options, 'Unoriented' );
else
    unoriented = false;
end

% Add lights to the axes
if isfield( Options, 'AddLights' )
    addLights = Options.AddLights;
    Options = rmfield( Options, 'AddLights' );
else
    addLights = false;
end

% Physical mesh vertex unit normals
if isfield( Options, 'VertexNormals' )
    VN = Options.VertexNormals;
    Options.rmfield( Options, 'VertexNormals');
    
    if ~isequal(size(VN), (VV))
        error('texture_patch_3d:inputs', ...
                'Invalid vertex normal input');
    end
else
    VN = per_vertex_normals( VV, FF, 'Weighting', 'angle' );
end
        

% Validate input texture image volume -------------------------------------
if (nargin < 5), IV = []; end

if ~isempty(IV) % Allow users to supply pre-made interpolant
    
    if isRGB

        if ~( iscell(IV) && (numel(IV) == 3) )
            
            error('texture_patch_3d:inputs', ...
                'Invalid texture image input');
            
        else
            
            goodIV = (ndims(IV{1}) == 3) && ...
                isequal(size(IV{1}), size(IV{2})) && ...
                isequal(size(IV{1}), size(IV{3})) && ...
                isequal(size(IV{2}), size(IV{3}));
            
            if ~goodIV
                error('texture_patch_3d:inputs', ...
                    'Invalid texture image input');
            end
            
        end
        
    else
        
        if ~(ndims(IV) == 3)
            error('texture_patch_3d:inputs', ...
            'Invalid texture image input');
        end
        
    end
    
end

% Create the interpolant object from the input image object ---------------
if isfield( Options, 'Interpolant' )
    
    IVIr = Options.Interpolant;
    
    if isRGB
        
        if ~( iscell(IVIr) && (numel(IVIr) == 3) )
            error('texture_patch_3d:inputs', ...
                'Invalid texture image interpolation object');
        end
        
        IVIg = IVIr{2}; IVIb = IVIr{3}; IVIr = IVIr{1};
        
    end
        
else
    
    if isRGB
        
        IVIr = griddedInterpolant(single(IV{1}), 'cubic');
        IVIg = griddedInterpolant(single(IV{2}), 'cubic');
        IVIb = griddedInterpolant(single(IV{3}), 'cubic');
        
    else
    
        IVIr = griddedInterpolant(single(IV), 'cubic');
        
    end
    
end

%--------------------------------------------------------------------------
% SCALAR FIELD COLORIZATION HANDLING
%--------------------------------------------------------------------------

if colorize
    
    % Map scalar field values to colors via the color mapping
    S(S < SCLim(1)) = SCLim(1);
    S(S > SCLim(2)) = SCLim(2);
    
    SInd = round( ( (S - SCLim(1)) ./ diff(SCLim) ) .* ...
        (size(SCMap,1)-1) ) + 1;
    SInd(isnan(SInd)) = 1;
    
    SColors = SCMap(SInd, :);
    
    % For colorizing grayscale images
    if ~isRGB, boneMap = bone (256); end
 
    
end

%--------------------------------------------------------------------------
% AMBIENT OCCLUSION HANDLING
%--------------------------------------------------------------------------

if applyAO
    
    if isempty(AO)
        
        % Calculate ambient occulsion
        % NOTE: MATLAB uses backwards normals for graphics...
        AO = ambient_occlusion(VV, FF, VV, -VN, AOSamples);
        if unoriented
            AO = min( AO, ambient_occlusion(VV, FF, VV, VN, AOSamples) );
        end
        
    end
    
end

%--------------------------------------------------------------------------
% FIND PATCH INTERPOLATION VALUES
%--------------------------------------------------------------------------

% The patch used for every triangle
% For simplicity these will be the same for each triangle regardless of
% size in either the physical or texture space. A smarter algorithm would
% probably take these features into account
Jr = zeros( (sizep+1), (sizep+1), 'single' );
if isRGB, Jg = Jr; Jb = Jr; end

% Linear indices of the 2D image associated with the triangle patch
jind = (sizep+1)^2:-1:1;

% Barycentric interpolation values for query points in the triangle patch
[ lambda1, lambda2, lambda3 ] = ...
    calculateBarycentricInterpolationValues( sizep );

%--------------------------------------------------------------------------
% CREATE THE SURFACE PLOT
%--------------------------------------------------------------------------

hold on

% Create container for all surface objects 
container = hggroup(gca) ;

% Loop through all of the faces of the mesh
for i = 1:size(FF,1)
    
    % The current face's physical vertex coordinates
    V = VV( FF(i,:), : );
    
    % The current face's texture vertex coordinates
    tV = TV( TF(i,:), : );
    
     % Define the physical face in 'surface' format
    x=[V(1,1) V(2,1); V(3,1) V(3,1)];
    y=[V(1,2) V(2,2); V(3,2) V(3,2)];
    z=[V(1,3) V(2,3); V(3,3) V(3,3)];
    
    % Define the texture coordinates of the 'surface'
    xyz = [ tV(1,:); tV(2,:); tV(3,:); tV(3,:) ];
    
    % Assemble vertex normal list
    vn = [ VN( FF(i,:), : ); VN( FF(i,3), : ) ];
    vn = cat(3, reshape(vn(:,1), [2 2]), reshape(vn(:,2), [2 2]), ...
        reshape(vn(:,3), [2 2]));
    
    % Calculate the texture interpolation coordinatex ---------------------
    
    pos(:,1) = xyz(1,1)*lambda1 + xyz(2,1)*lambda2 + xyz(3,1)*lambda3;
    pos(:,2) = xyz(1,2)*lambda1 + xyz(2,2)*lambda2 + xyz(3,2)*lambda3;
    pos(:,3) = xyz(1,3)*lambda1 + xyz(2,3)*lambda2 + xyz(3,3)*lambda3;
    
    % Map texture to surface image ----------------------------------------
    Jr(jind) = IVIr( pos(:,1), pos(:,2), pos(:,3) ) ;
    if Imin > -Inf
        Jr(Jr < Imin) = Imin ;
    end
    if Imax < Inf 
        Jr(Jr > Imax) = Imax ;
        Jr = Jr / Imax ;
    end
    J(:,:,1) = Jr;
    if isRGB
        Jg(jind) = IVIg( pos(:,1), pos(:,2), pos(:,3) );
        Jb(jind) = IVIb( pos(:,1), pos(:,3), pos(:,3) );
        if Imin > -Inf
            Jg(Jg < Imin) = Imin ;
            Jb(Jb < Imin) = Imin ;
        end
        if Imax < Inf 
            Jg(Jg > Imax) = Imax ;
            Jb(Jb > Imax) = Imax ;
            Jb = Jb / Imax ;
            Jg = Jg / Imax ;
        end
        J(:,:,2) = Jg;
        J(:,:,3) = Jb;
    end
    
    % Apply affine transformations ----------------------------------------
    if rotate
        xyz = [x(1) y(1) z(1); x(3) y(3) z(3); x(2) y(2) z(2)] ;
        xyzp = (rot * xyz')' ;
        x = [xyzp(1, 1) xyzp(2, 1); xyzp(3, 1) xyzp(3, 1)];
        y = [xyzp(1, 2) xyzp(2, 2); xyzp(3, 2) xyzp(3, 2)];
        z = [xyzp(1, 3) xyzp(2, 3); xyzp(3, 3) xyzp(3, 3)];
    end
    
    if translate
        x = x + trans(1) ;
        y = y + trans(2) ;
        z = z + trans(3) ;
    end
    
    if dilate
        x = x * dilation ;
        y = y * dilation ;
        z = z * dilation ;
    end
    
    % Apply scalar field colorization -------------------------------------
    
    if colorize
        
        % Grayscale texture mappings must be made RGB 
        if ~isRGB
            
            cind = round( 255 .* (Jr(:)-min(Jr(:))) ./ range(Jr(:)) ) + 1;
            cind(isnan(cind)) = 1;
            
            Jr = reshape(boneMap(cind, 1), sizep+1, sizep+1);
            Jg = reshape(boneMap(cind, 2), sizep+1, sizep+1);
            Jb = reshape(boneMap(cind, 3), sizep+1, sizep+1);
            
        end
        
        if numel(S) == size(FF,1) % Face based scalar field
            
            SC = SColors(i,:); % Face has a single flat color
            
            % Calculate alpha blending for each channel
            Jr = SAlpha .* SC(1) + (1-SAlpha) .* Jr;
            Jg = SAlpha .* SC(2) + (1-SAlpha) .* Jg;
            Jb = SAlpha .* SC(3) + (1-SAlpha) .* Jb;
            
        else % Vertex based scalar field
            
            % Interpolated scalar field colors from vertices
            SCr = SColors(FF(i,1),1)*lambda1 ...
                + SColors(FF(i,2),1)*lambda2 + SColors(FF(i,3),1)*lambda3;
            SCg = SColors(FF(i,1),2)*lambda1 ...
                + SColors(FF(i,2),2)*lambda2 + SColors(FF(i,3),2)*lambda3;
            SCb = SColors(FF(i,1),3)*lambda1 ...
                + SColors(FF(i,2),3)*lambda2 + SColors(FF(i,3),3)*lambda3;
            
            % Reshape to the size of the surface image
            SCr = reshape(SCr, sizep+1, sizep+1);
            SCg = reshape(SCg, sizep+1, sizep+1);
            SCb = reshape(SCb, sizep+1, sizep+1);
            
            % Calculate alpha blending for each channel
            Jr = SAlpha .* SCr + (1-SAlpha) .* Jr;
            Jg = SAlpha .* SCg + (1-SAlpha) .* Jg;
            Jb = SAlpha .* SCb + (1-SAlpha) .* Jb;
            
        end
        
        % Combine channels
        J = cat(3, Jr, Jg, Jb);
        
    end
    
    % Apply ambient occlusion ---------------------------------------------
    if applyAO
        
        % Interpolate the ambient occlusion values on vertices
        AOMat = AO(FF(i,1))*lambda1 + AO(FF(i,2))*lambda2 + ...
            AO(FF(i,3))*lambda3;
        
        % Reshape to the size of the surface image
        AOMat = reshape(AOMat, sizep+1, sizep+1);
        
        Jr = (1-AOFactor) .* Jr + AOFactor .* (Jr .* (1-AOMat));
        J(:,:,1) = Jr;
        if (isRGB || colorize)
            Jg = (1-AOFactor) .* Jg + AOFactor .* (Jg .* (1-AOMat));
            Jb = (1-AOFactor) .* Jb + AOFactor .* (Jb .* (1-AOMat));
            J(:,:,2) = Jg;
            J(:,:,3) = Jb;
        end
        
    end
        
    % Show surface --------------------------------------------------------
    % NOTE: MATLAB uses left-handed normals for graphics objects...
    Options.VertexNormals = -vn;
    % J 
    surface( container, x, y, z, J, Options );
    
end

% t1 = hgtransform('Parent',gca);
% Txy = makehgtform('translate',Options.translate);
% set(t1,'Matrix',Txy)

hold off

%--------------------------------------------------------------------------
% Lighting Handling
%--------------------------------------------------------------------------

if addLights
    
    % Extract the center of the current axis bounding box
    XLim = get(gca, 'XLim');
    YLim = get(gca, 'YLim');
    ZLim = get(gca, 'ZLim');
    
    cen = [ mean(XLim) mean(YLim) mean(ZLim) ];
    
    light( 'Position', [cen(1) cen(2) 10*(max(ZLim)-min(ZLim))+cen(3)], ...
        'Style', 'local', 'Color', [1 1 1]/3 );
    light( 'Position', [cen(1) 10*(max(ZLim)-min(ZLim))+cen(2) cen(3)], ...
        'Style', 'local', 'Color', [1 1 1]/3 );
    light( 'Position', [cen(1) 10*(min(ZLim)-max(ZLim))+cen(2) cen(3)], ...
        'Style', 'local', 'Color', [1 1 1]/3 );
    
end

end


%==========================================================================
% CALCULATEBARYCENTRICINTERPOLATIONVALUES
%==========================================================================

function [ lambda1, lambda2, lambda3 ] = ...
    calculateBarycentricInterpolationValues( sizep )
% CALCULATEBARYCENTRICINTERPOLATIONVALUES Calculates the barycentric
% interpolation values of a set of query points in a square of size
% (sizep+1)X(sizep+1) relative to an abstract triangle in the upper region
% of that square

% Define the triangle in the upperpart of the square
x1 = sizep; y1 = sizep;
x2 = sizep; y2 = 0;
x3 = 0; y3 = sizep;

% Calculate the barycentric coordinates for each query point on the square
% (These will correspond to the linear indices corresponding to the 2D
% image associated with the triangle)
[x,y] = ndgrid(0:sizep,0:sizep); % NOTE THE USE OF NDGRID
x = x(:); y = y(:);

detT = (x1-x3)*(y2-y3) - (x2-x3)*(y1-y3);
lambda1 = ( (y2-y3).*(x-x3) + (x3-x2).*(y-y3) ) ./ detT;
lambda2 = ( (y3-y1).*(x-x3) + (x1-x3).*(y-y3) ) ./ detT;
lambda3 = 1 - lambda1 - lambda2;

end


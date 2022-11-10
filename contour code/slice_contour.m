function contour_slice = slice_contour(contour,dimension,value)

% calculates vertices of a slice through a contour

NRP=size(contour,1);
contour_slice=cell(NRP,1);

for ind_RP=1:NRP
    verts_in=contour{ind_RP};

    if ~isempty(verts_in)
        % get list of vertices of simplices describing contour surface
        ndim=size(verts_in,2);
        TRI=convhulln(verts_in);
        verts=zeros(length(TRI),ndim,ndim);
        for i=1:ndim
            verts(:,:,i)=verts_in(TRI(:,i),:);
        end

        % restrict to only simplices that contain target value
        Xmin = min(squeeze(verts(:,dimension,:)),[],2);
        Xmax = max(squeeze(verts(:,dimension,:)),[],2);
        good = Xmin<=value & Xmax>=value;
        verts = verts(good,:,:);

        % interpolate edges of simplices
        ntri=length(verts);
        perms=nchoosek(1:ndim,2);
        nperm=size(perms,1);
        verts_interp = nan(ntri,ndim-1,nperm);
        dims=1:ndim;
        dims(dimension)=[];
        for i=1:nperm
            v1=squeeze(verts(:,:,perms(i,1)));
            v2=squeeze(verts(:,:,perms(i,2)));
            frac=(v2(:,dimension)-value)./(v2(:,dimension)-v1(:,dimension));
            bad=frac<0 | frac>1;
            frac(bad)=NaN;
            ind=0;
            for j=dims
                ind=ind+1;
                verts_interp(:,ind,i) = frac.*v1(:,j) + (1-frac).*v2(:,j);
            end
        end

        % keep only unique vertices
        verts_slice = zeros(ntri*nperm,ndim-1);
        for i=1:ndim-1
            v=squeeze(verts_interp(:,i,:));
            verts_slice(:,i)=v(:);
        end
        bad=isnan(verts_slice(:,1));
        verts_slice(bad,:)=[];
        verts_slice=unique(verts_slice,'rows');

        % set outputs
        contour_slice{ind_RP}=verts_slice;
    end
end

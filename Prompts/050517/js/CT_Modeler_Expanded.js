
/////// CHALKTALK MODELER FOR HTML5 ///////////////////////////////////////////////////////////////

var isPhone = function(){ return false; };

if (! window.CT) CT = { };
CT.REVISION = "0";

CT.def = function(a, b) { return a !== undefined ? a : b !== undefined ? b : 0; }
CT.fovToFL = function(fov) { return 1 / Math.tan(fov / 2); }
CT.time = 0;
CT.imu = { compass:null, alpha:null, beta:null, gamma:null, ax:null, ay:null, az:null };
// if (window.DeviceOrientationEvent)
//    window.addEventListener('deviceorientation', function(event) {
//       CT.imu.alpha = event.alpha;  // COMPASS DIRECTION IN DEGREES
//       CT.imu.beta  = event.beta;   // TILT FRONT-2-BACK IN DEGREES
//       CT.imu.gamma = event.gamma;  // TILT LEFT-2-RIGHT IN DEGREES

//       if (CT.imu.compass = event.webkitCompassHeading) // IF COMPASS HEADING IS SUPPORTED,
//          CT.imu.alpha = CT.imu.compass;                // THEN USE IT.

//       if (CT.imu.callback)
//          CT.imu.callback();
//    });
// if (window.DeviceMotionEvent)
//    window.addEventListener('devicemotion', function(event) {
//       CT.imu.ax = event.acceleration.x;
//       CT.imu.ay = event.acceleration.y;
//       CT.imu.az = event.acceleration.z;
//    });


CT.cross     = function(a, b) { return [ a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0] ]; }
CT.dot       = function(a, b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
CT.normalize = function(v)    { var d = Math.sqrt(CT.dot(v,v)); v[0]/=d; v[1]/=d; v[2]/=d; return v; }

CT.copy      = function(m, src)   { for(var i = 0 ; i<m.length ; i++) m[i] = src[i];  return m; }
CT.identity  = function(m)        { for(var i = 0 ; i<16 ; i++) m[i] = i % 5 ? 0 : 1; return m; }
CT.rotateX   = function(m, a)     { return CT.matrixMultiply(m, CT.matrixRotatedX  (a)    , m); }
CT.rotateY   = function(m, a)     { return CT.matrixMultiply(m, CT.matrixRotatedY  (a)    , m); }
CT.rotateZ   = function(m, a)     { return CT.matrixMultiply(m, CT.matrixRotatedZ  (a)    , m); }
CT.scale     = function(m, x,y,z) { return CT.matrixMultiply(m, CT.matrixScaled    (x,y,z), m); }
CT.translate = function(m, x,y,z) { return CT.matrixMultiply(m, CT.matrixTranslated(x,y,z), m); }

CT.matrixToQ = function(m) {
   if (m[0] + m[5] + m[10] > 0) {        var s = 2 * Math.sqrt(m[0] + m[5] + m[10] + 1);
                                         return [ (m[9] - m[6]) / s, (m[2] - m[8]) / s, (m[4] - m[1]) / s, s / 4 ];
   } else if (m[0]>m[5] && m[0]>m[10]) { var s = 2 * Math.sqrt(1 + m[0] - m[5] - m[10]);
                                         return [ s / 4, (m[1] + m[4]) / s, (m[2] + m[8]) / s, (m[9] - m[6]) / s ];
   } else if (m[5]>m[10]) {              var s = 2 * Math.sqrt(1 + m[5] - m[0] - m[10]);
                                         return [ (m[1] + m[4]) / s, s / 4, (m[6] + m[9]) / s, (m[2] - m[8]) / s ];
   } else {                              var s = 2 * Math.sqrt(1 + m[10] - m[0] - m[5]);
                                         return [ (m[2] + m[8]) / s, (m[6] + m[9]) / s, s / 4, (m[4] - m[1]) / s ];
} }
CT.matrixFromPQ = function(pq)      { var qx = pq[3], qy = pq[4], qz = pq[5], qw = pq[6];
                                      return [ 1 - 2*qy*qy - 2*qz*qz,     2*qx*qy - 2*qz*qw,     2*qz*qx + 2*qy*qw, 0,
                                                   2*qx*qy + 2*qz*qw, 1 - 2*qx*qx - 2*qz*qz,     2*qy*qz - 2*qx*qw, 0,
                                                   2*qz*qx - 2*qy*qw,     2*qy*qz + 2*qx*qw, 1 - 2*qx*qx - 2*qy*qy, 0,
                                               pq[0], pq[1], pq[2], 1 ]; }
CT.matrixIdentity = function()      { return [1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1]; }
CT.matrixRotatedX = function(a)     { return [1,0,0,0, 0,Math.cos(a),Math.sin(a),0, 0,-Math.sin(a),Math.cos(a),0, 0,0,0,1]; }
CT.matrixRotatedY = function(a)     { return [Math.cos(a),0,-Math.sin(a),0, 0,1,0,0, Math.sin(a),0,Math.cos(a),0, 0,0,0,1]; }
CT.matrixRotatedZ = function(a)     { return [Math.cos(a),Math.sin(a),0,0, -Math.sin(a),Math.cos(a),0,0, 0,0,1,0, 0,0,0,1]; }
CT.matrixScaled   = function(x,y,z) { if (x instanceof Array) { z = x[2]; y = x[1]; x = x[0]; }
                                      return [x,0,0,0, 0,CT.def(y,x),0,0, 0,0,CT.def(z,x),0, 0,0,0,1]; }
CT.matrixTranslated = function(x,y,z) { if (x instanceof Array) { z = x[2]; y = x[1]; x = x[0]; }
                                        return [1,0,0,0, 0,1,0,0, 0,0,1,0, x,y,z,1]; }

CT.matrixMultiply = function(a, b, dst) {
   var tmp = [], n;
   for (n = 0 ; n < 16 ; n++)
      tmp.push(a[n&3] * b[n&12] + a[n&3 | 4] * b[1 | n&12] + a[n&3 | 8] * b[2 | n&12] + a[n&3 | 12] * b[3 | n&12]);
   if (dst)
      CT.copy(dst, tmp);
   return CT.def(dst, tmp);
}
CT.matrixInverse = function(src) {
  function s(col, row) { return src[col & 3 | (row & 3) << 2]; }
  function cofactor(c0, r0) {
     var c1 = c0+1, c2 = c0+2, c3 = c0+3, r1 = r0+1, r2 = r0+2, r3 = r0+3;
     return (c0 + r0 & 1 ? -1 : 1) * ( (s(c1, r1) * (s(c2, r2) * s(c3, r3) - s(c3, r2) * s(c2, r3)))
                                     - (s(c2, r1) * (s(c1, r2) * s(c3, r3) - s(c3, r2) * s(c1, r3)))
                                     + (s(c3, r1) * (s(c1, r2) * s(c2, r3) - s(c2, r2) * s(c1, r3))) );
  }
  var n, dst = [], det = 0;
  for (n = 0 ; n < 16 ; n++) dst.push(cofactor(n >> 2, n & 3));
  for (n = 0 ; n <  4 ; n++) det += src[n] * dst[n << 2];
  for (n = 0 ; n < 16 ; n++) dst[n] /= det;
  return dst;
}
CT.matrixTransform = function(m, v)  {
    var x = v[0], y = v[1], z = v[2], w = CT.def(v[3], 1);
    return [ x*m[0] + y*m[4] + z*m[ 8] + w*m[12],
             x*m[1] + y*m[5] + z*m[ 9] + w*m[13],
             x*m[2] + y*m[6] + z*m[10] + w*m[14],
             x*m[3] + y*m[7] + z*m[11] + w*m[15] ]; }
CT.matrixTranspose = function(m) { return [m[0],m[4],m[ 8],m[12],m[1],m[5],m[ 9],m[13],
                                           m[2],m[6],m[10],m[14],m[3],m[7],m[11],m[15]]; }

CT.boxesVertices = function(B) {
   var H, L, i, v = [];
   for (i = 0 ; i < B.length ; i++) {
      L = B[i][0];
      H = B[i][1];
      v.push( L[0], L[1], L[2],   H[0], L[1], L[2],   L[0], H[1], L[2],   H[0], H[1], L[2],
              L[0], L[1], H[2],   H[0], L[1], H[2],   L[0], H[1], H[2],   H[0], H[1], H[2] );
   }
   return v;
}
CT.boxesFaces = function(B) {
   return CT.shapesFaces(B, [ 0,2,3, 3,1,0,   4,5,7, 7,6,4,   0,4,6, 6,2,0,
                              1,3,7, 7,5,1,   0,1,5, 5,4,0,   2,6,7, 7,3,2 ]);
}
CT.shapesFaces = function(S, V) {
   var f = [], i, n, nv = 0;
   for (n = 0 ; n < V.length ; n++)
      nv = max(nv, V[n] + 1);
   for (i = 0 ; i < S.length ; i++)
      for (n = 0 ; n <= V.length ; n++)
         f.push(nv * i + V[n]);
   return f;
}

CT.Scene = function(canvas) {
   this._fog = [0,0,0,0];
   this._fov = Math.PI / 3;
   this._iod = 0.0625;
   this._gl = canvas.getContext('experimental-webgl');
   this._lColor = [0,0,0, 0,0,0, 0,0,0];
   this._lDirInfo = [0,0,0, 0,0,0, 0,0,0];
   this._objects = [];
   this._viewMatrix = [1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,this.getFL(),1];
}
CT.Scene.prototype = {
   add      : function(obj)   { this._objects.push(obj); return obj._scene = this; },
   getFL    : function()      { return CT.fovToFL(this._fov); },
   getFOV   : function()      { return this._fov; },
   getIOD   : function()      { return this._iod; },
   getHDir  : function()      { this._updateLightVectors(); return this._hDir; },
   getLDir  : function()      { this._updateLightVectors(); return this._lDir; },
   getObj   : function(i)     { return this._objects[i]; },
   getStereo: function()      { return CT.def(this._stereo); },
   setFog   : function(fog)   { this._fog = fog; return this; },
   setFOV   : function(fov)   { this._fov = fov; return this; },
   setIOD   : function(iod)   { this._iod = iod; return this; },
   setStereo: function(y_n)   { this._stereo = y_n ? 1 : 0; return this; },
   doDepthTest : function(yes) { // TOGGLE WHETHER TO DO A DEPTH TEST WHEN DRAWING
      var gl = this._gl;
      if (yes) {
         gl.enable(gl.DEPTH_TEST);
         gl.depthFunc(gl.LEQUAL);
      }
      else
         gl.disable(gl.DEPTH_TEST);
   },
   getViewMatrix : function() {
      return this._viewMatrix;
   },
   getViewMatrixInverse : function() {
      if (! this._viewMatrixInverse)
         this._viewMatrixInverse = CT.matrixInverse(this._viewMatrix);
      return this._viewMatrixInverse;
   },
   remove : function(obj) {
      for (var i = 0 ; i < this._objects.length ; i++)
         if (obj == this._objects[i]) {
	    this._objects.splice(i, 1);
	    break;
	 }
   },
   setCanvas : function(canvas) {
      this._gl = canvas.getContext('experimental-webgl');
      for (var i = 0 ; i < this._objects.length ; i++)
         this._objects[i].setGL(this._gl);
   },
   setLight : function(n, lDir, lColor) {
      CT.normalize(lDir);
      for (var i = 0 ; i < 3 ; i++) {
	 this._lDirInfo[3 * n + i] = lDir[i];
         this._lColor  [3 * n + i] = lColor ? lColor[i] : 1;
      }
      delete this._lDir;
      return this;
   },
   setViewMatrix : function(matrix) {
      CT.copy(this._viewMatrix, matrix);
      delete this._viewMatrixInverse;
      delete this._lDir;
      return this;
   },
   _updateLightVectors : function() {
      if (! this._lDir) {
         var m = this.getViewMatrixInverse(), v = this._lDirInfo, dir, i;
         this._lDir = [];
         this._hDir = [];
         for (i = 0 ; i < v.length ; i += 3) {
            dir = CT.normalize(CT.matrixTransform(m, [ v[i], v[i+1], v[i+2], 0 ]));
	    this._lDir.push(dir[0], dir[1], dir[2]);
            dir = CT.normalize([dir[0], dir[1], dir[2] + 1]);
	    this._hDir.push(dir[0], dir[1], dir[2]);
         }
      }
   },
}

CT.Object = function() { } // BASE CLASS FOR COMPOSITE 3D OBJECTS.

CT.Object.prototype = {
   clear            : function()        { this._children = []; return this; },
   getChild         : function(i)       { return this._children[i]; },
   getProperty      : function(p, dflt) { return this[p] ? this[p] : this._parent ? this._parent.getProperty(p,dflt) : dflt; },
   getScene         : function()        { return this._scene = this.getProperty('_scene'); },
   getShape         : function()        { return this._shape; },
   identity         : function()        { this._matrix = CT.matrixIdentity(); return this; },
   numChildren      : function()        { return this._children.length; },
   rotateX          : function(a)       { return this.transform(CT.matrixRotatedX(a)); },
   rotateY          : function(a)       { return this.transform(CT.matrixRotatedY(a)); },
   rotateZ          : function(a)       { return this.transform(CT.matrixRotatedZ(a)); },
   scale            : function(x, y, z) { return this.transform(CT.matrixScaled(x, y, z)); },
   setColor         : function(r,g,b,p) { if (r instanceof Array) { p = CT.def(r[3]); b = r[2]; g = r[1]; r = r[0]; }
                                              return this.setPhong([.1*r,.1*g,.1*b, .5*r,.5*g,.5*b, .5,.5,.5, CT.def(p,20), 1]); },
   setFragmentShader: function(str)     { this._fragmentShader = str; return this; },
   setMatrix        : function(src)     { CT.copy(this._matrix, src); return this; },
   setMetal         : function(r,g,b,p) { return this.setPhong([r/80,g/80,b/80, r/20,g/20,b/20, r,g,b, CT.def(p,7), 1]); },
   setNormalMap     : function(file)    { this._textureFile1 = file; return this; },
   setOpacity       : function(t)       { this._opacity = t; return this; },
   setPhong         : function(phong)   { if (! this._phong) ; this._phong = phong; return this; },
   setPQ            : function(pq)      { return this.setMatrix(CT.matrixFromPQ(pq)); },
   setTexture       : function(file)    { this._textureFile0 = file; return this; },
   setVertexShader  : function(str)     { this._vertexShader = str; return this; },
   translate        : function(x, y, z) { return this.transform(CT.matrixTranslated(x,y,z)); },
   transform        : function(m)       { CT.matrixMultiply(this._matrix, m, this._matrix); return this; },

   addChild : function(child) {
      this._children.push(child);
      child._parent = this;
      return child;
   },
   draw : function(globalMatrix) {
      if (! globalMatrix) {
	 var scene = this.getScene();
         if (this._gl != scene._gl)
	    this.setGL(scene._gl);
         CT.copy(this._viewMatrixInverse, scene.getViewMatrixInverse());
	 globalMatrix = this._viewMatrixInverse;
	 if (CT.imu.alpha != null) {
	    if (CT.imu.alpha0 === undefined) {
	       CT.imu.alpha0 = CT.imu.alpha;
	       this.getScene().setStereo(true);
            }

            if (! CT.imu.matrix)
	       CT.imu.matrix = CT.matrixIdentity();

	    CT.identity(CT.imu.matrix);
            if (window.tracked)
               CT.translate(CT.imu.matrix, -.1 * tracked[1][0], -.1 * tracked[1][1], -.1 * (tracked[1][2] - 2 * scene.getFL()));
	    CT.rotateX  (CT.imu.matrix,                  CT.imu.gamma  * Math.PI / 180 + Math.PI/2);
	    CT.rotateZ  (CT.imu.matrix,                  CT.imu.beta   * Math.PI / 180            );
	    CT.rotateY  (CT.imu.matrix, (CT.imu.alpha0 - CT.imu.alpha) * Math.PI / 180            );

	    CT.matrixMultiply(globalMatrix, CT.imu.matrix, globalMatrix);

            window.imuData = CT.imu;

            if (window.server)
	       server.broadcastGlobal('imuData');
         }
      }
      CT.matrixMultiply(globalMatrix, this._matrix, this._globalMatrix);
      this._drawShape(this._shape, this._globalMatrix);
      for (var i = 0 ; i < this._children.length ; i++) 
          this._children[i].draw(this._globalMatrix);
   },
   getFocus : function() {
      var m = this.getGlobalMatrix();
      return -m[14] / Math.sqrt(m[12] * m[12] + m[13] * m[13] + m[14] * m[14]);
   },
   getGlobalMatrix : function() {
      return this._globalMatrix ? this._globalMatrix : CT.matrixIdentity();
   },
   init : function(shape) {
      this._children = [];
      this.identity();
      this._viewMatrixInverse = CT.matrixIdentity();
      this._globalMatrix = CT.matrixIdentity();
      this._shape = shape;
      if (shape)
         shape._object = this;
      return this;
   },
   setGL : function(gl) {
      this._gl = gl;
      if (this._shape && this._shape._gl != gl)
         this._shape._initGL(gl);
      this._loadTexture(this._textureFile0, '_texture0');
      this._loadTexture(this._textureFile1, '_texture1');
      for (var i = 0 ; i < this._children.length ; i++)
         this._children[i].setGL(gl);
   },
   _drawShape : function(shape, globalMatrix) {
      if (shape) {
         shape._renderMatrix = globalMatrix;
         shape._useProgram();
         shape._draw();
      }
   },
   _loadTexture : function(fileName, textureId) {
      if (fileName && ! this[textureId]) {
         var that = this, gl = this._gl, image = new Image(), texture = gl.createTexture();
         image.onload = function() {
	    try {
               gl.activeTexture (gl.TEXTURE0);
               gl.bindTexture   (gl.TEXTURE_2D, texture);
               gl.pixelStorei   (gl.UNPACK_FLIP_Y_WEBGL, true);
               gl.texImage2D    (gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, image);
               gl.texParameteri (gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
               gl.texParameteri (gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR_MIPMAP_NEAREST);
               gl.generateMipmap(gl.TEXTURE_2D);
               that[textureId] = texture;
            }
	    catch(e) { console.log(e); }
         };
         image.src = fileName;
      }
   },
}

CT.Shape = function() { } // A SHAPE IS THE OBJECT THAT HANDLES RENDERING OF A SINGLE TRIANGLE STRIP.

CT._vertexShaderPerspective =
   ['   pos.y *= uAspect;'
   ,'   pos = mix(pos, vec4(.5 * pos.xyz, pos.w), uEye * uEye);'
   ,'   pos.x += pos.w * uEye * .5;'
   ,'   vClipX = uEye * pos.x;'
   ,'   pos.z = pos.z * .1;'
   ,'   gl_Position = pos;'
   ].join('\n');

CT.Shape.prototype = {
   getProperty     : function(p,dflt) { return this._object.getProperty(p, dflt); },
   getScene        : function()       { return this._object.getScene(); },
   getSceneProperty: function(p)      { return (this.getScene())[p]; },
   _address        : function(name)   { return this._gl.getUniformLocation(this._program, name); },

   surfaceExtruded : function(nu, nv, fu, fv) {
      this.surfaceParametric(nu, nv, function(u, v) {
         var xy = fu(u, v), x = xy[0], y = xy[1], p = fv(v), p1 = fv(v+.001);
	 var dz = CT.normalize([ p1[0]-p[0], p1[1]-p[1], p1[2]-p[2] ]);
	 var dx = [1, 0, 0];
	 var dy = CT.normalize(CT.cross(dz, dx));
         return [ p[0] + x*dx[0] + y*dy[0]  ,  p[1] + x*dx[1] + y*dy[1]  ,  p[2] + x*dx[2] + y*dy[2] ];
      });
   },
   surfaceRevolved : function(nu, nv, f) {
      this.surfaceParametric(nu, nv, function(u, v) {
         var theta = 2 * Math.PI * u, xz = f(v);
         return [ xz[0] * Math.cos(theta), xz[0] * Math.sin(theta), xz[1] ];
      });
   },
   surfaceParametric : function(nu, nv, f) {              // CONVERT USER FUNCTION TO A PARAMETRIC SURFACE.
      var du = 1 / nu, dv = 1 / nv;
      var vertices = [], s;
      var addVertex = function(u, v) {                                             // COMPUTE PARAMETRIC VERTEX.
         var p = f( Math.max(0, Math.min(1, u)), Math.max(0, Math.min(.999, v)) ), // f MUST EVAL TO A POINT.
             pu = f(u + du/200, v), nu = f(u - du/200, v),                         // APPROXIMATE TANGENTS BY
             pv = f(u, v + dv/200), nv = f(u, v - du/200);                         // FINITE DIFFERENCING.
         var uu = CT.normalize([ pu[0] - nu[0], pu[1] - nu[1], pu[2] - nu[2] ]);
         var vv = CT.normalize([ pv[0] - nv[0], pv[1] - nv[1], pv[2] - nv[2] ]);
         var ww = CT.cross(uu, vv);
         vertices.push(p[0],p[1],p[2],  uu[0],uu[1],uu[2],  vv[0],vv[1],vv[2],  ww[0],ww[1],ww[2],  u,v,  0,0);
      }
      var u = 0, d = du;
      for (var v = 0 ; v < 1 ; v += dv, d = -d)               // ZIGZAG ACROSS ROWS TO FORM A TRIANGLE STRIP.
         for (var i = 0 ; i <= nu ; i++) {
            addVertex(u, v);
	    addVertex(u, v + dv);
	    if (i < nu)
	       u += d;
            else
               addVertex(u, v);
         }
      this._vertices = new Float32Array(vertices);
   },
   _addShader : function(type, str) {
      var gl = this._gl;
      var s = gl.createShader(type);
      gl.shaderSource(s, str);
      gl.compileShader(s);
      if (! gl.getShaderParameter(s, gl.COMPILE_STATUS))
        console.log("shader compiler error: " + gl.getShaderInfoLog(s));
      gl.attachShader(this._program, s);
   },
   _attrib : function(name, size, loc, stride) {
      var gl = this._gl, bpe = Float32Array.BYTES_PER_ELEMENT;
      gl.enableVertexAttribArray(this[name]);
      gl.vertexAttribPointer(this[name], size, gl.FLOAT, false, stride * bpe,  loc * bpe);
   },
   _cameraMatrix : function(f, eye) {
      return [ 1,0,0,0,  0,1,0,0,  0,0,0,-1/f,  -.5*eye,0,-1,0 ];
   },
   _draw : function() {
      var gl = this._gl, program = this._program, scene = this.getScene();
      gl.uniform4fv (this._address('uFog'    ), scene._fog);
      gl.uniform3fv (this._address('uHDir'   ), scene.getHDir());
      gl.uniform3fv (this._address('uLColor' ), this._lColor);
      gl.uniform3fv (this._address('uLDir'   ), scene.getLDir());
      gl.uniform1f  (this._address('uOpacity'), this.getProperty('_opacity', 1));
      gl.uniform1fv (this._address('uPhong'  ), this.getProperty('_phong', [.05,.05,.05, .5,.5,.5, 1,1,1, 20, 1]));
      gl.uniform1fv (this._address('uTexture'), [this.getProperty('_texture0') ? 1 : 0 ,
                                                 this.getProperty('_texture1') ? 1 : 0 ] );
      gl.uniformMatrix4fv(this._address('uNorMatrix'), false,
         (function(m) { var x = m[0] * m[0] + m[1] * m[1] + m[ 2] * m[ 2],   // THE NORMAL MATRIX
                            y = m[4] * m[4] + m[5] * m[5] + m[ 6] * m[ 6],   // IS THE INVERSE
                            z = m[8] * m[8] + m[9] * m[9] + m[10] * m[10];   // TRANSPOSE OF THE
                        return [ m[0]/x, m[1]/x, m[ 2]/x, 0,                 // POSITION MATRIX.
                                 m[4]/y, m[5]/y, m[ 6]/y, 0,
                                 m[8]/z, m[9]/z, m[10]/z, 0,  0,0,0,1 ]; })(this._renderMatrix));

      this._drawArrays(this._vertices.length / 16);
   },
   _drawArrays : function(len) {
      var eye, fl, gl, matrix, stereo;
      fl = this.getScene().getFL();
      iod = this.getScene().getIOD();
      gl = this._gl;
      var drawMode = this instanceof CT.ShapePolyhedron ? gl.TRIANGLES : gl.TRIANGLE_STRIP;
      stereo = this.getScene().getStereo();
      gl.uniform1f(this._address('uAspect'), gl.canvas.width / gl.canvas.height);
      gl.uniform1f(this._address('uTime'), CT.time);
      for (eye = -stereo ; eye <= stereo ; eye += 2) {
         gl.uniformMatrix4fv(this._address('uMatrix'), false,
            CT.matrixMultiply(this._cameraMatrix(fl, eye * iod), this._renderMatrix));
         gl.uniform1f(this._address('uEye'), eye);
	 gl.drawArrays(drawMode, 0, len);
      }
   },
   _initGL : function(gl) {
      this._gl           = gl;
      this._lColor       = this.getSceneProperty('_lColor');
      this._lDir         = this.getSceneProperty('_lDir');
      this._hDir         = this.getSceneProperty('_hDir');
      this._renderMatrix = CT.matrixIdentity();

      gl.enable    (gl.DEPTH_TEST);
      gl.depthFunc (gl.LEQUAL);
      gl.enable    (gl.BLEND);
      gl.blendFunc (gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);

      this._program = gl.createProgram();
      this._addShader(gl.VERTEX_SHADER  , this.getProperty('_vertexShader'  , this._defaultVertexShader  ));
      this._addShader(gl.FRAGMENT_SHADER, this.getProperty('_fragmentShader', this._defaultFragmentShader));
      gl.linkProgram(this._program);
      this._buffer  = gl.createBuffer();

      this._aPos      = gl.getAttribLocation(this._program, 'aPos');
      this._aTangent  = gl.getAttribLocation(this._program, 'aTangent');
      this._aBinormal = gl.getAttribLocation(this._program, 'aBinormal');
      this._aNormal   = gl.getAttribLocation(this._program, 'aNormal');
      this._aUV       = gl.getAttribLocation(this._program, 'aUV');
   },
   _isLoaded : function() {
      var gl = this._gl, texture, i;
      if (gl._program == this._program)
         return true;
      gl._program = this._program;
      gl.useProgram(this._program);

      gl.bindBuffer(gl.ARRAY_BUFFER, this._buffer);
      gl.bufferData(gl.ARRAY_BUFFER, this._vertices, gl.STATIC_DRAW);
      gl.uniform1iv(this._address('uSampler'), [0,1,2,3,4,5,6,7]);
      for (i = 0 ; i < 8 ; i++)
         if (texture = this.getProperty('_texture' + i)) {
            gl.activeTexture(gl.TEXTURE0 + i);
            gl.bindTexture  (gl.TEXTURE_2D, texture);
         }
      return false;
   },
   _useProgram : function() {
      if (! this._isLoaded()) {
         this._attrib('_aPos'     , 3,  0, 16);
         this._attrib('_aTangent' , 3,  3, 16);
         this._attrib('_aBinormal', 3,  6, 16);
         this._attrib('_aNormal'  , 3,  9, 16);
         this._attrib('_aUV'      , 2, 12, 16);
      }
   },
   _defaultVertexShader :
      ['precision highp float;'
      ,'attribute vec3 aPos, aTangent, aBinormal, aNormal;'
      ,'attribute vec2 aUV;'
      ,'uniform mat4   uNorMatrix, uMatrix;'
      ,'uniform float  uEye, uAspect, uTime;'
      ,'varying vec4   vPos;'
      ,'varying vec3   vTangent, vBinormal, vNormal;'
      ,'varying vec2   vUV;'
      ,'varying float  vClipX;'
      ,'void main() {'
      ,'   vTangent  = normalize(uNorMatrix * vec4(aTangent ,0.)).xyz;'
      ,'   vBinormal = normalize(uNorMatrix * vec4(aBinormal,0.)).xyz;'
      ,'   vNormal   = normalize(uNorMatrix * vec4(aNormal  ,0.)).xyz;'
      ,'   vUV = aUV;'
      ,'   vec4 pos = uMatrix * vec4(aPos, 1.);'
      ,'   vPos = pos;'
      ,    CT._vertexShaderPerspective
      ,'}'
      ].join('\n'),

   _defaultFragmentShader :
      ['precision highp float;'
      ,'uniform sampler2D uSampler[8];'
      ,'uniform vec4  uFog;'
      ,'uniform vec3  uLColor[3], uLDir[3], uHDir[3];'
      ,'uniform float uEye, uOpacity, uPhong[10], uTexture[2];'
      ,'varying vec4  vPos;'
      ,'varying vec3  vTangent, vBinormal, vNormal;'
      ,'varying vec2  vUV;'
      ,'varying float vClipX;'
      ,'void main(void) {'
      ,'   if (vClipX < 0.) discard;'
      ,'   vec4 b = texture2D(uSampler[1], vUV);'
      ,'   vec3 normal = normalize(mix(vNormal, b.r * vTangent + b.g * vBinormal + b.b * vNormal, b.a * uTexture[1]));'
      ,'   vec3 rgb = vec3(uPhong[0], uPhong[1], uPhong[2]);'
      ,'   for (int i = 0 ; i < 3 ; i++) {'
      ,'      float d = max(0., dot(uLDir[i], normal));'
      ,'      float s = max(0., dot(uHDir[i], normal));'
      ,'      rgb += uLColor[i] * (vec3(uPhong[3], uPhong[4], uPhong[5]) * d +'
      ,'                           vec3(uPhong[6], uPhong[7], uPhong[8]) * pow(s, 4. * uPhong[9]));'
      ,'   }'
      ,'   vec4 texture = texture2D(uSampler[0], vUV);'
      ,'   rgb = mix(rgb, rgb * texture.rgb, texture.a * uTexture[0]);'
      ,'   rgb = mix(uFog.rgb, rgb, exp(-uFog.a * (2. * vPos.w - 1.)));'
      ,'   gl_FragColor = vec4(sqrt(rgb), uOpacity);'
      ,'}'
      ].join('\n'),
}

CT.ShapePolyhedron = function(vs, fs, ns, parity) {
   this._vertices = (function() {
      function addTriangle(a, b, c, p) {
         var ax = vs[3 * a], ay = vs[3 * a + 1], az = vs[3 * a + 2],
             bx = vs[3 * b], by = vs[3 * b + 1], bz = vs[3 * b + 2],
             cx = vs[3 * c], cy = vs[3 * c + 1], cz = vs[3 * c + 2],
             U = CT.normalize(            [bx - ax, by - ay, bz - az]),
             W = CT.normalize(CT.cross(U, [cx - bx, cy - by, cz - bz])),
             V = CT.normalize(CT.cross(W, U));

         if (p === undefined)
	    p = 0;
         var q = 1 - p;

         if (ns) {
            vertices.push(ax,ay,az, U[0],U[1],U[2], V[0],V[1],V[2], ns[3 * a],ns[3 * a + 1],ns[3 * a + 2], p,p, 0,0);
            vertices.push(bx,by,bz, U[0],U[1],U[2], V[0],V[1],V[2], ns[3 * b],ns[3 * b + 1],ns[3 * b + 2], p,q, 0,0);
            vertices.push(cx,cy,cz, U[0],U[1],U[2], V[0],V[1],V[2], ns[3 * c],ns[3 * c + 1],ns[3 * c + 2], q,q, 0,0);
         }
         else {
            vertices.push(ax,ay,az, U[0],U[1],U[2], V[0],V[1],V[2], W[0],W[1],W[2], p,p, 0,0);
            vertices.push(bx,by,bz, U[0],U[1],U[2], V[0],V[1],V[2], W[0],W[1],W[2], p,q, 0,0);
            vertices.push(cx,cy,cz, U[0],U[1],U[2], V[0],V[1],V[2], W[0],W[1],W[2], q,q, 0,0);
         }
      }
      var vertices = [], parity = 0;
      for (var i = 0 ; i < fs.length ; i += 3)
            addTriangle(fs[i], fs[i + 1], fs[i + 2], parity = 1 - parity);
      return new Float32Array(vertices);
   })();
}
CT.ShapePolyhedron.prototype = new CT.Shape;

CT.ShapeCube = function() {
   this._vertices = (function() {
      function addFace(x,y,z, ax,ay,az, bx,by,bz) {
         vertices.push(x-ax-bx, y-ay-by, z-az-bz,  ax,ay,az, bx,by,bz,  x,y,z,  0,0, 0,0);
         vertices.push(x+ax-bx, y+ay-by, z+az-bz,  ax,ay,az, bx,by,bz,  x,y,z,  1,0, 0,0);
         vertices.push(x+ax+bx, y+ay+by, z+az+bz,  ax,ay,az, bx,by,bz,  x,y,z,  1,1, 0,0);
         vertices.push(x-ax+bx, y-ay+by, z-az+bz,  ax,ay,az, bx,by,bz,  x,y,z,  0,1, 0,0);
         vertices.push(x-ax-bx, y-ay-by, z-az-bz,  ax,ay,az, bx,by,bz,  x,y,z,  0,0, 0,0);
      }
      var vertices = [];
      addFace( 1, 0, 0,   0, 1, 0,   0, 0, 1);
      addFace( 0, 1, 0,   0, 0, 1,   1, 0, 0);
      addFace( 0, 0, 1,   1, 0, 0,   0, 1, 0);
      addFace( 0, 0,-1,  -1, 0, 0,   0, 1, 0);
      addFace( 0,-1, 0,   0, 0,-1,   1, 0, 0);
      addFace(-1, 0, 0,   0,-1, 0,   0, 0, 1);
      return new Float32Array(vertices);
   })();
}

CT.LINE_BUF_SIZE = 200;

CT.ShapePath = function() {
   this.draw = function() {
      this.setRGBA([1,0,0,1]);
      this.prototype.draw();
   }
   this._draw = function() {
      var path = this.getProperty('_path');
      for (var i = 4 ; i < path.length ; i += 4) // REMOVE ANY REPEATED POINTS BEFORE RENDERING.
         if ( path[i    ] == path[i - 4] &&
              path[i + 1] == path[i - 3] &&
              path[i + 2] == path[i - 2] &&
        path[i + 3] == path[i - 1] ) {
            path.splice(i, 4);
      i -= 4;
   }
      if (! path)
         return;
      var gl      = this._gl,
          isFill  = this.getProperty('_isFill', false),
          program = this._program,
          rgba    = this.getProperty('_rgba', [1,1,1,1]);

      if (window.displayListener)
         this._object.sendToServer(isFill, path, rgba);
         //server.broadcastObject( { display:'path', isFill: isFill, path: path, rgba: rgba } );

      gl.uniform1f       (this._address('uFill'  ), isFill ? 1 : 0);
      gl.uniform4fv      (this._address('uRgba'     ), [pow(rgba[0],.45),
                                                        pow(rgba[1],.45),
                                                        pow(rgba[2],.45), rgba[3]]);

      if (isFill) {                                   //------------- FILL ---------------
         var centroid = [0,0,0];
         for (var i = 0 ; i < path.length ; i += 4)
             for (var j = 0 ; j < 3 ; j++)
                centroid[j] += path[i + j] / (path.length / 4);
         gl.uniform3fv(this._address('uCentroid'), centroid);

   var p, p0 = path.slice(0, 4), p1;            // REMEMBER VERY FIRST POINT.
         for (var i = 0 ; path.length ; i++) {        // KEEP GOING UNTIL PATH IS EMPTY.
      p = path.splice(0, CT.LINE_BUF_SIZE);     // EACH SEGMENT IS A BUFFERFUL OF POINTS.
      if (i > 0)
         p = p1.concat(p);                      // PREPEND LAST POINT FROM PREV SEGMENT, IF ANY.
      p1 = p.slice(p.length - 4, p.length);     // REMEMBER LAST POINT OF THIS SEGMENT.
      if (path.length == 0)
         p = p.concat(p0);                      // APPEND VERY FIRST POINT TO LAST SEGMENT.
            gl.uniform4fv(this._address('uPath'), p);
            this._drawArrays(2 * (p.length / 4));
         }
      }
      else {                                          //------------- STROKE --------------
         while (path.length) {
      var p = path.splice(0, CT.LINE_BUF_SIZE);
            gl.uniform1f (this._address('uNpts'   ), p.length / 4);
            gl.uniform4fv(this._address('uPath'   ), p);
            gl.uniform1f (this._address('uRadius0'), this.getProperty('_radius0', 1));
            gl.uniform1f (this._address('uRadius1'), this.getProperty('_radius1', 1));
            this._drawArrays(2 * (p.length / 4 + 10));
      if (path.length > 0) {
         var n = p.length;
         path.unshift(p[n-4], p[n-3], p[n-2], p[n-1]);
            }
         }
      }
   };
   this._useProgram = function() {
      if (! this._isLoaded())
   this._attrib('_aPos', 2, 0);
   };
   this._vertices = (function() {
      var vertices = [];
      for (var y = 0 ; y < CT.LINE_BUF_SIZE ; y++)
         for (var x = 0 ; x <= 1 ; x++)
            vertices.push(x, y);
      return new Float32Array(vertices);
   })();
   this._defaultVertexShader =
      ['precision highp float;'
      ,'attribute vec2 aPos;'
      ,'uniform mat4  uMatrix;'
      ,'uniform vec4  uPath[' + CT.LINE_BUF_SIZE + '];'
      ,'uniform vec3  uCentroid;'
      ,'uniform float uFill, uNpts, uRadius0, uRadius1, uEye, uAspect, uTime;'
      ,'varying float vClipX;'
      ,'vec4 P(int i) { return uMatrix * vec4(uPath[i].xyz, 1.); }'
      ,'void main() {'
      ,'   vec4 pos;'
      ,'   if (uFill > 0.5)'
      ,'      pos = uMatrix * vec4(mix(uCentroid, uPath[int(aPos.y)].xyz, aPos.x), 1.);'
      ,'   else {'
      ,'      float x = mix(-1., 1., aPos.x);'
      ,'      float radius = mix(uRadius0, uRadius1, clamp((aPos.y - 5.) / uNpts, 0., 1.));'
      ,'      vec2 dp;'
      ,'      int n;'
      ,'      if (aPos.y <= 5.) {'
      ,'         n = 0;'
      ,'         pos = P(n);'
      ,'         vec4 p1 = mix(pos, P(n+1), 0.01);'
      ,'         vec2 u = normalize(vec2(pos.y - p1.y, p1.x - pos.x));'
      ,'         vec2 v = normalize(vec2(p1.x - pos.x, p1.y - pos.y));'
      ,'         float a = 3.14159 / 2. * aPos.y / 5.;'
      ,'         dp = x * sin(a) * u - cos(a) * v;'
      ,'      }'
      ,'      else if (aPos.y >= uNpts + 5.) {'
      ,'         n = int(uNpts) - 1;'
      ,'         pos = P(n);'
      ,'         vec4 p1 = mix(pos, P(n-1), -0.01);'
      ,'         vec2 u = normalize(vec2(pos.y - p1.y, p1.x - pos.x));'
      ,'         vec2 v = normalize(vec2(p1.x - pos.x, p1.y - pos.y));'
      ,'         float a = 3.14159 / 2. * (uNpts + 10. - aPos.y) / 5.;'
      ,'         dp = x * sin(a) * u + cos(a) * v;'
      ,'      }'
      ,'      else {'
      ,'         n = int(min(uNpts - 2., aPos.y - 5.));'
      ,'         pos = P(n);'
      ,'         vec4 p1 = mix(pos, P(n+1), 0.01);'
      ,'         float dx = p1.x - pos.x;'
      ,'         float dy = p1.y - pos.y;'
      ,'         dx = sign(dx) * max(abs(dx), 0.00001);'
      ,'         dy = sign(dy) * max(abs(dy), 0.00001);'
      ,'         dp = x * normalize(vec2(-dy, dx));'
      ,'      }'
      ,'      pos += vec4(radius * dp * uPath[n].w, 0., 0.);'
      ,'   }'
      ,    CT._vertexShaderPerspective
      ,'}'
      ].join('\n');

   this._defaultFragmentShader =
      ['precision highp float;'
      ,'uniform vec4  uRgba;'
      ,'varying float vClipX;'
      ,'void main() {'
      ,'   if (vClipX < 0.) discard;'
      ,'   gl_FragColor = uRgba;'
      ,'}'
     ].join('\n');
}; CT.ShapePath.prototype = new CT.Shape;

CT.Path = function() {

// TODO: COMBINE MULTIPLE CALLS TO SINGLE draw() CALLS WHEREVER POSSIBLE
//       THAT IS: WHEN rgba MATCHES AND NOT isFill.

   this.drawsPerFrame = 0;

   this.drawPath = function() {
      this.draw();
      this.drawsPerFrame++;
   };

   this.sendToServer = function(isFill, path, rgba) {
      if (this._pathsData === undefined)
         this._pathsData = [];

      this._pathsData.push([isFill, path, rgba]);
   };

   this.flush = function() {
/*
      if (window.displayListener) {
         server.broadcastObject( { display: 'flush' } );
*/
      var dataArray, i, n, npts, pd;

      if (window.displayListener) {
         pd = this._pathsData;

   npts = 0;
   for (var n = 0 ; n < pd.length ; n++)
      npts += pd[n][1].length;

         dataArray = [];
   dataArray.push(0);
   dataArray.push(pd.length);

   for (n = 0 ; n < pd.length ; n++) {
      dataArray.push(pd[n][0] << 24 | pd[n][1].length);
      dataArray.push(pd[n][2]);
   }

   for (n = 0 ; n < pd.length ; n++)
      for (i = 0 ; i < pd[n][1].length ; i++)
         dataArray.push(pd[n][1][j]);

   dataArray[0] = dataArray.length;

         server.broadcastObject( { display: dataArray } );
         this._pathsData = [];
      }

      this.drawsPerFrame = 0;
   };

   this.setFill = function(isFill) {
      this._isFill = isFill;
   };
   this.setLineWidth = function(width0, width1) {
      var scale = isPhone() ? 1 : 0.5;
      this._radius0 = scale * width0;
      this._radius1 = scale * def(width1, width0);
   };
   this.setPath = function(src) {
      function P(_i) { return path[path.length - _i]; }
      function setP(_i, value) { path[path.length - _i] = value; }

      function addPointsAtSharpBends(src) {
         path = [];
         for (var i = 0 ; i < src.length ; i += 4) {
            var x = src[i], y = src[i + 1], z = src[i + 2], w = src[i + 3];
            if (path.length) {
               var ax = P(4) - P(8), ay = P(3) - P(7), az = P(2) - P(6),
             bx = x    - P(4), by = y    - P(3), bz = z    - P(2);
         var aa = ax * ax + ay * ay + az * az,
             ab = ax * bx + ay * by + az * bz,
       bb = bx * bx + by * by + bz * bz;
         if (ab * ab < .9 * aa * bb) {
                  path.push(P(4), P(3), P(2), P(1));
                  setP(8, P(8) - .01 * ax);                // INSERT EXTRA POINTS
                  setP(7, P(7) - .01 * ay);                // AT SHARP BENDS
                  setP(6, P(6) - .01 * az);                // IN THE PATH.
      if (i < src.length - 4) {
                     x += .01 * bx;
                     y += .01 * by;
                     z += .01 * bz;
                  }
         }
            }
            path.push(x, y, z, w);
         }
   return path;
      }

      var path = this._isFill ? src : addPointsAtSharpBends(src);
      if (path.length == 8) {
         path.push(P(4), P(3), P(2), P(1));
         path.push(P(4), P(3), P(2), P(1));
         for (var j = 1 ; j <= 4 ; j++) {
            setP(8 + j, mix(P(12 + j), P(j), .01));
            setP(4 + j, mix(P(12 + j), P(j), .99));
         }
      }
      this._path = path;
   };
   this.setRGBA = function(rgba) {
      this._rgba = cloneArray(rgba);
   };
   this.init(new CT.ShapePath());
};
CT.Path.prototype = new CT.Object;

CT.ShapeCube.prototype = new CT.Shape;

(CT.ShapeDisk         = function(n) { this.surfaceRevolved(n,1,function(t) { return [t+.001,0]; }); }).prototype = new CT.Shape;
(CT.ShapeExtruded     = function(nu,nv,fu,fv) { this.surfaceExtruded(nu, nv, fu, fv); }).prototype = new CT.Shape;
(CT.ShapeOpenCylinder = function(n) { this.surfaceRevolved(n, 2, function(t) { return [1, 2*t-1] }); }).prototype = new CT.Shape;
(CT.ShapeParametric   = function(nu,nv,f) { this.surfaceParametric(nu, nv, f); }).prototype = new CT.Shape;
(CT.ShapeRevolved     = function(nu,nv,f) { this.surfaceRevolved(nu, nv, f); }).prototype = new CT.Shape;

(CT.ShapeSphere = function(nu, nv, s) {
   nu = nu ? nu : 24;
   nv = nv ? nv : Math.floor(nu/2);
   s = CT.def(s);
   this.surfaceRevolved(nu, nv, function(t) { var phi = Math.PI*(Math.max(s, t) - .5); return [ Math.cos(phi), Math.sin(phi) ]; });
}).prototype = new CT.Shape;

(CT.ShapeSquare = function(nu, nv) {
   this.surfaceParametric(1, 1, function(u, v) { return [2*u-1, 2*v-1, 0]; });
}).prototype = new CT.Shape;

(CT.ShapeTorus = function(nu, nv, r) {
   if (! r) r = 0.3;
   this.surfaceRevolved(nu, nv, function(t) { var phi = 2*Math.PI*t; return [ 1 - r * Math.cos(phi), -r * Math.sin(phi) ]; });
}).prototype = new CT.Shape;

// HERE ARE THE TYPES OF OBJECTS THAT THE MODELER CURRENTLY SUPPORTS.

CT.Cube         = function()        { this.init(new CT.ShapeCube());               }; CT.Cube.prototype         = new CT.Object;
CT.Cylinder     = function(n)       { this.init();
                                      if (n === undefined) n = 16;
                                      this.addChild(new CT.OpenCylinder(n));
                                      this.addChild(new CT.Disk(n)).translate(0,0,-1);
                                      this.addChild(new CT.Disk(n)).translate(0,0,1);
                                                                                   }; CT.Cylinder.prototype     = new CT.Object;
CT.Disk         = function(n)       { this.init(new CT.ShapeDisk(n));              }; CT.Disk.prototype         = new CT.Object;
CT.Extruded = function(nu,nv,fu,fv) { this.init(new CT.ShapeExtruded(nu,nv,fu,fv));}; CT.Extruded.prototype     = new CT.Object;
CT.Node         = function()        { this.init();                                 }; CT.Node.prototype         = new CT.Object;
CT.OpenCylinder = function(n)       { if (n === undefined) n = 16;
                                      this.init(new CT.ShapeOpenCylinder(n));      }; CT.OpenCylinder.prototype = new CT.Object;
CT.Parametric   = function(nu,nv,f) { this.init(new CT.ShapeParametric(nu,nv,f));  }; CT.Parametric.prototype   = new CT.Object;
CT.Polyhedron   = function(vs,fs,ns){ this.init(new CT.ShapePolyhedron(vs,fs,ns)); }; CT.Polyhedron.prototype   = new CT.Object;
CT.Revolved     = function(nu,nv,f) { this.init(new CT.ShapeRevolved(nu,nv,f));    }; CT.Revolved.prototype     = new CT.Object;
CT.Sphere       = function(m,n,s)   { this.init(new CT.ShapeSphere(m,n,s));        }; CT.Sphere.prototype       = new CT.Object;
CT.Square       = function()        { this.init(new CT.ShapeSquare());             }; CT.Square.prototype       = new CT.Object;
CT.Torus        = function(m,n,r)   { this.init(new CT.ShapeTorus(m,n,r));         }; CT.Torus.prototype        = new CT.Object;

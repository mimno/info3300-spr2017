var Modeler = new (function(){
	var scene, obj;
	var index = 0;
	var rotationSpeed = [ 0, 0, 0 ];


	this.init = function(elementId){
		scene = new CT.Scene(document.getElementById(elementId));
        scene.setLight(0, [1,1,1]);

        obj = new CT.Node();
        scene.add(obj);

	    var lastTime = new Date().getTime() / 1000;
	    setTimeout(function() {
	        function update() {
		        for (var i = 0 ; i < obj.numChildren() ; i++){
		            var child = obj.getChild(i);
		            // child.identity();
		            // .translate(4*(i%4)-6, i<4?2:-2, 0).rotateY(time).rotateX(time/2);
		        }

		        obj.draw();
		    }

	        setInterval(function() {
	            window.time = (new Date()).getTime() / 1000;
	            var deltaTime = window.time - lastTime;
	            lastTime = window.time;

	            obj.rotateX(deltaTime * rotationSpeed[0]);
	            obj.rotateY(deltaTime * rotationSpeed[1]);
	            obj.rotateZ(deltaTime * rotationSpeed[2]);
	            update();

	            ++index;
	        }, 16);
	    }, 100);
	}

	this.addPoint = function(x, y, z, color){
		var sphere = new CT.Sphere(16, 8);
        sphere.translate(x, y, z);
        sphere.scale(0.01, 0.01, 0.01);
        color = color || [ 1, 0, 0 ];
        color.push(0.5);
        sphere.setColor(color);
        obj.addChild(sphere);
	}

	this.scale = function(scalar){
		obj.scale(scalar);
	}

	this.setRotationSpeed = function(array){
		while(array.length < 3){
			array.push(0);
		}
		rotationSpeed = array;
	}

	this.setRotation = function(array){
		while(array.length < 3){
			array.push(0);
		}
		obj.identity().rotateX(array[0]).rotateY(array[1]).rotateZ(array[2]);
	}
})();
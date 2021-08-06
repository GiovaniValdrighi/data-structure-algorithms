




function update(){
    var urls = {'cat': {'random': {
        '100': "https://raw.githubusercontent.com/GiovaniValdrighi/data-structure-algorithms/main/delaunay-triangulation/images/cat_100_random_points.json",
        "1000": "https://raw.githubusercontent.com/GiovaniValdrighi/data-structure-algorithms/main/delaunay-triangulation/images/cat_1000_random_points.json",
        "10000": "https://raw.githubusercontent.com/GiovaniValdrighi/data-structure-algorithms/main/delaunay-triangulation/images/cat_10000_random_points.json",
        }, "edge" :  {
            '100': "https://raw.githubusercontent.com/GiovaniValdrighi/data-structure-algorithms/main/delaunay-triangulation/images/cat_100_points.json",
            "1000": "https://raw.githubusercontent.com/GiovaniValdrighi/data-structure-algorithms/main/delaunay-triangulation/images/cat_1000_points.json",
            "10000": "https://raw.githubusercontent.com/GiovaniValdrighi/data-structure-algorithms/main/delaunay-triangulation/images/cat_10000_points.json"
        }
    },
    'face': {'random': {
        '100': "https://raw.githubusercontent.com/GiovaniValdrighi/data-structure-algorithms/main/delaunay-triangulation/images/giovani_100_random_points.json",
        "1000": "https://raw.githubusercontent.com/GiovaniValdrighi/data-structure-algorithms/main/delaunay-triangulation/images/giovani_1000_random_points.json",
        "10000": "https://raw.githubusercontent.com/GiovaniValdrighi/data-structure-algorithms/main/delaunay-triangulation/images/giovani_10000_random_points.json",
        }, "edge" :  {
            '100': "https://raw.githubusercontent.com/GiovaniValdrighi/data-structure-algorithms/main/delaunay-triangulation/images/giovani_100_points.json",
            "1000": "https://raw.githubusercontent.com/GiovaniValdrighi/data-structure-algorithms/main/delaunay-triangulation/images/giovani_1000_points.json",
            "10000": "https://raw.githubusercontent.com/GiovaniValdrighi/data-structure-algorithms/main/delaunay-triangulation/images/giovani_10000_points.json"
        }
    }
}

    var rad1 = document.f1.picture
    var rad2 = document.f2.points
    var rad3 = document.f3.n_points
    var rad4 = document.f4.color
    //var rad5 = document.f5.triangulation
    var picture, points, n_points, color//, triangulation
    for(let i = 0; i < rad1.length; i++){
        if(rad1[i].checked){
            picture = rad1[i].value
        }
    }
    for(let i = 0; i < rad2.length; i++){
        if(rad2[i].checked){
            points = rad2[i].value
        }
    }
    for(let i = 0; i < rad3.length; i++){
        if(rad3[i].checked){
            n_points = rad3[i].value
        }
    }
    for(let i = 0; i < rad4.length; i++){
        if(rad4[i].checked){
            color = rad4[i].value
        }
    }
    //for(let i = 0; i < rad5.length; i++){
    //    if(rad5[i].checked){
    //        triangulation = rad5[i].value;
    //    }
    //}
   
    var url = urls[picture][points][n_points]
    fetch(url)
    .then(function(response){
        return response.json()
    }).then(function(response){
        console.log(response)
        var colors = response.map(d => d.colors)

        //normalize coordinates
        var coords = [].concat(...response.map(d => d.points))
        var coords_norm = []
        var xmax = Math.max(...coords.map(d => d[0]))
        var xmin = Math.min(...coords.map(d => d[0]))
        var ymax = Math.max(...coords.map(d => d[1]))
        var ymin = Math.min(...coords.map(d => d[1]))
        for(let i = 0; i < coords.length; i++){
            coords_norm.push([(coords[i][0] - xmin)/(xmax - xmin)*2 - 1,
                            1 - (coords[i][1] - ymin)/(ymax - ymin)*2 ])
        }

        //normalize colors
        var colors_norm = []
        var colors = [].concat(...response.map(d => d.colors))
        if(color == "mean"){
            for(let i = 0; i < (colors.length/3); i++){
                //var b = false
                //if(triangulation == 'holes'){
                //    b = Math.random() < 0.1
                //    console.log(b)
                //}
                for(let j = 0; j < 3; j++){
                    var mean = (colors[3*i][j] + colors[3*i + 1][j] + colors[3*i+2][j])/3 
                    //mean = b ? 239 : mean
                    colors[3*i][j] = mean
                    colors[3*i + 1][j] = mean
                    colors[3*i + 2][j] = mean
                }
            }
        }
        
        for(let i = 0; i < colors.length; i++){
            var c = colors[i].map(d => d/255).concat([1])
            colors_norm.push(c);
        }
        console.log(colors_norm)

        document.getElementById("plot").innerHTML = "";
        var canvas = document.createElement('canvas');
        canvas.width = parseInt(xmax - xmin)
        canvas.height = parseInt(ymax - ymin)

        document.getElementById("plot").appendChild(canvas);
        var regl = createREGL(canvas)

        function drawTriangle(n){
            return  regl({
                vert : `
                precision mediump float;
                attribute vec2 position;
                attribute vec4 color;
                varying vec4 fragColor;
                void main () {
                gl_Position = vec4(position, 0, 1);
                fragColor = color;
                }`,
        
                frag : `
                precision mediump float;
                varying vec4 fragColor;
                void main(){
                gl_FragColor = fragColor; //vec4(1, 0, 0, 1);//
                }
                `,
            
                attributes : {
                position: regl.prop("position"),
                color: regl.prop("color")
                },
            
                count: n,
                primitive: 'triangle'
            })
        } 

        drawTriangle(3*response.length)({
            position: coords_norm,
            color: colors_norm
        })
        
    })

}


update()

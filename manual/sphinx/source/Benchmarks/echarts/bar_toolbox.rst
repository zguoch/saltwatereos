.. raw:: html

	<body>
	    <div id="0fb67cd45e2c4a1499ab87543d21e5a3" class="chart-container" style="width:900px; height:500px;"></div>
	    <script>
	        var chart_0fb67cd45e2c4a1499ab87543d21e5a3 = echarts.init(
	            document.getElementById('0fb67cd45e2c4a1499ab87543d21e5a3'), 'white', {renderer: 'canvas'});
	        var option_0fb67cd45e2c4a1499ab87543d21e5a3 = {
	    "animation": true,
	    "animationThreshold": 2000,
	    "animationDuration": 1000,
	    "animationEasing": "cubicOut",
	    "animationDelay": 0,
	    "animationDurationUpdate": 300,
	    "animationEasingUpdate": "cubicOut",
	    "animationDelayUpdate": 0,
	    "color": [
	        "#c23531",
	        "#2f4554",
	        "#61a0a8",
	        "#d48265",
	        "#749f83",
	        "#ca8622",
	        "#bda29a",
	        "#6e7074",
	        "#546570",
	        "#c4ccd3",
	        "#f05b72",
	        "#ef5b9c",
	        "#f47920",
	        "#905a3d",
	        "#fab27b",
	        "#2a5caa",
	        "#444693",
	        "#726930",
	        "#b2d235",
	        "#6d8346",
	        "#ac6767",
	        "#1d953f",
	        "#6950a1",
	        "#918597"
	    ],
	    "series": [
	        {
	            "type": "bar",
	            "name": "\u5546\u5bb6A",
	            "legendHoverLink": true,
	            "data": [
	                93,
	                116,
	                47,
	                107,
	                58,
	                91,
	                141
	            ],
	            "showBackground": false,
	            "barMinHeight": 0,
	            "barCategoryGap": "20%",
	            "barGap": "30%",
	            "large": false,
	            "largeThreshold": 400,
	            "seriesLayoutBy": "column",
	            "datasetIndex": 0,
	            "clip": true,
	            "zlevel": 0,
	            "z": 2,
	            "label": {
	                "show": true,
	                "position": "top",
	                "margin": 8
	            }
	        },
	        {
	            "type": "bar",
	            "name": "\u5546\u5bb6B",
	            "legendHoverLink": true,
	            "data": [
	                20,
	                77,
	                64,
	                71,
	                35,
	                32,
	                56
	            ],
	            "showBackground": false,
	            "barMinHeight": 0,
	            "barCategoryGap": "20%",
	            "barGap": "30%",
	            "large": false,
	            "largeThreshold": 400,
	            "seriesLayoutBy": "column",
	            "datasetIndex": 0,
	            "clip": true,
	            "zlevel": 0,
	            "z": 2,
	            "label": {
	                "show": true,
	                "position": "top",
	                "margin": 8
	            }
	        }
	    ],
	    "legend": [
	        {
	            "data": [
	                "\u5546\u5bb6A",
	                "\u5546\u5bb6B"
	            ],
	            "selected": {
	                "\u5546\u5bb6A": true,
	                "\u5546\u5bb6B": true
	            },
	            "show": false,
	            "padding": 5,
	            "itemGap": 10,
	            "itemWidth": 25,
	            "itemHeight": 14
	        }
	    ],
	    "tooltip": {
	        "show": true,
	        "trigger": "item",
	        "triggerOn": "mousemove|click",
	        "axisPointer": {
	            "type": "line"
	        },
	        "showContent": true,
	        "alwaysShowContent": false,
	        "showDelay": 0,
	        "hideDelay": 100,
	        "textStyle": {
	            "fontSize": 14
	        },
	        "borderWidth": 0,
	        "padding": 5
	    },
	    "xAxis": [
	        {
	            "show": true,
	            "scale": false,
	            "nameLocation": "end",
	            "nameGap": 15,
	            "gridIndex": 0,
	            "inverse": false,
	            "offset": 0,
	            "splitNumber": 5,
	            "minInterval": 0,
	            "splitLine": {
	                "show": false,
	                "lineStyle": {
	                    "show": true,
	                    "width": 1,
	                    "opacity": 1,
	                    "curveness": 0,
	                    "type": "solid"
	                }
	            },
	            "data": [
	                "\u5c0f\u7c73",
	                "\u4e09\u661f",
	                "\u534e\u4e3a",
	                "\u82f9\u679c",
	                "\u9b45\u65cf",
	                "VIVO",
	                "OPPO"
	            ]
	        }
	    ],
	    "yAxis": [
	        {
	            "show": true,
	            "scale": false,
	            "nameLocation": "end",
	            "nameGap": 15,
	            "gridIndex": 0,
	            "inverse": false,
	            "offset": 0,
	            "splitNumber": 5,
	            "minInterval": 0,
	            "splitLine": {
	                "show": false,
	                "lineStyle": {
	                    "show": true,
	                    "width": 1,
	                    "opacity": 1,
	                    "curveness": 0,
	                    "type": "solid"
	                }
	            }
	        }
	    ],
	    "title": [
	        {
	            "text": "Bar-\u663e\u793a ToolBox",
	            "padding": 5,
	            "itemGap": 10
	        }
	    ],
	    "toolbox": {
	        "show": true,
	        "orient": "horizontal",
	        "itemSize": 15,
	        "itemGap": 10,
	        "left": "80%",
	        "feature": {
	            "saveAsImage": {
	                "type": "png",
	                "backgroundColor": "auto",
	                "connectedBackgroundColor": "#fff",
	                "show": true,
	                "title": "\u4fdd\u5b58\u4e3a\u56fe\u7247",
	                "pixelRatio": 1
	            },
	            "restore": {
	                "show": true,
	                "title": "\u8fd8\u539f"
	            },
	            "dataView": {
	                "show": true,
	                "title": "\u6570\u636e\u89c6\u56fe",
	                "readOnly": false,
	                "lang": [
	                    "\u6570\u636e\u89c6\u56fe",
	                    "\u5173\u95ed",
	                    "\u5237\u65b0"
	                ],
	                "backgroundColor": "#fff",
	                "textareaColor": "#fff",
	                "textareaBorderColor": "#333",
	                "textColor": "#000",
	                "buttonColor": "#c23531",
	                "buttonTextColor": "#fff"
	            },
	            "dataZoom": {
	                "show": true,
	                "title": {
	                    "zoom": "\u533a\u57df\u7f29\u653e",
	                    "back": "\u533a\u57df\u7f29\u653e\u8fd8\u539f"
	                },
	                "icon": {},
	                "xAxisIndex": false,
	                "yAxisIndex": false,
	                "filterMode": "filter"
	            },
	            "magicType": {
	                "show": true,
	                "type": [
	                    "line",
	                    "bar",
	                    "stack",
	                    "tiled"
	                ],
	                "title": {
	                    "line": "\u5207\u6362\u4e3a\u6298\u7ebf\u56fe",
	                    "bar": "\u5207\u6362\u4e3a\u67f1\u72b6\u56fe",
	                    "stack": "\u5207\u6362\u4e3a\u5806\u53e0",
	                    "tiled": "\u5207\u6362\u4e3a\u5e73\u94fa"
	                },
	                "icon": {}
	            },
	            "brush": {
	                "icon": {},
	                "title": {
	                    "rect": "\u77e9\u5f62\u9009\u62e9",
	                    "polygon": "\u5708\u9009",
	                    "lineX": "\u6a2a\u5411\u9009\u62e9",
	                    "lineY": "\u7eb5\u5411\u9009\u62e9",
	                    "keep": "\u4fdd\u6301\u9009\u62e9",
	                    "clear": "\u6e05\u9664\u9009\u62e9"
	                }
	            }
	        }
	    }
	};
	        chart_0fb67cd45e2c4a1499ab87543d21e5a3.setOption(option_0fb67cd45e2c4a1499ab87543d21e5a3);
	    </script>
	</body>

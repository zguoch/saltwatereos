// import the extension via require
var example = require("./example");

// // calling the global method
// var x = 42;
// var y = 105;
// var g = example.gcd(x, y);

// // Accessing the global variable
// var f = example.Foo;
// example.Foo = 3.1415926;
var aaa = example.fact(8);
var foo = example.Foo;
var f = new foo();
var c = f.bar(20);
console.log(aaa, c);

/**
 * 
 */
const http = require('http');

function dispatch(op, request, response) {
    switch (op) {
        case 'worker':
            console.log('worker');
            response['message'] = 'Worker';
            break;

        default:
            response['message'] = 'NO OP';
    }
}

http.createServer(function (req, res) {
    let op = req.url;
    if (op.endsWith("/"))
        op = op.substr(0, op.length - 1);
    let n = op.lastIndexOf("/");
    if (!(n < 0))
        op = op.substr(n + 1);

    let data = '';
    req.on('data', chunk => {
        data += chunk;
    })
    req.on('end', () => {
        let request = JSON.parse(data);
        let response = {};
        response['status'] = 'OK';

        console.log(request);
        dispatch(op, request, response);
        console.log(response);

        res.writeHead(200, { 'Content-Type': 'application/json', 'Access-Control-Allow-Origin': "*" });
        res.end(JSON.stringify(response));
    })
    req.on('error', err => {
        res.writeHead(500, { 'Content-Type': 'text/html' });
        res.end(error.message);
    })
}).listen(8080);

console.log("server started");

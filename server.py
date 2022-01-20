import socket
import struct
import numpy as np
import numpy.linalg

HOST = '127.0.0.1'
PORT = 8000

server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
server.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
server.bind((HOST, PORT))
server.listen(10)

while True:
    conn, addr = server.accept()
    # clientMessage = str(conn.recv(1024), encoding='utf-8')
    try:
        size = conn.recv(4)
        size, = struct.unpack('I', size)
        decode = ''
        for i in range(int(size*size)):
            decode = decode + 'd'
        print('size : ' + str(size*size))

        matrix = conn.recv(5000)
        print('Get matrix')
        data = list()
        data = struct.unpack(decode, matrix)
    
        data = np.array(data).reshape(size, size)

        # data = np.array(data, dtype=float)
        temp=numpy.linalg.pinv(data)
        pseudo_inverse = list()
        for row in temp:
            for element in row:
                pseudo_inverse.append(element)

        reply_matrix = struct.pack(decode, *pseudo_inverse)
        print('send back pinv!!\n')
        conn.sendall(reply_matrix)
        conn.close()
    except Exception as e:
        print(e)
        conn.close()
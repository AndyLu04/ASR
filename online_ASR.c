#include "lib.h"
#include "online_ASR.h"
#include "queue.h"
#include "ASR.h"
#include "filter.h"

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

void online_ASR()
{
    queue* my_queue = create_queue();

    pthread_t tcp_thread;
    pthread_t process_thread;

    pthread_create(&tcp_thread, NULL, TCP_client, my_queue);
    pthread_create(&process_thread, NULL, MW_ASR, my_queue);

    pthread_join(tcp_thread, NULL);
    pthread_join(process_thread, NULL);
}

void* TCP_client(void* argv)
{
    queue* my_queue = (queue*)argv;
    int sockfd = 0;
    sockfd = socket(AF_INET , SOCK_STREAM , 0);

    if (sockfd == -1){
        printf("Fail to create a socket.");
    }

    //socket connection
    struct sockaddr_in info;
    bzero(&info,sizeof(info));
    info.sin_family = PF_INET;

    //localhost
    info.sin_addr.s_addr = inet_addr("10.211.55.2");
    info.sin_port = htons(7000);

    int err = connect(sockfd,(struct sockaddr *)&info,sizeof(info));
    if(err==-1){
        printf("Connection error");
    }

    char server_reply[5000] = {};

    while(1)
    {
        recv(sockfd,server_reply,sizeof(server_reply),0);
        double* data = malloc(8 * sizeof(double));
        memcpy(data, server_reply, 8 * sizeof(double));
        pthread_mutex_lock(&mutex);
        push_queue(my_queue, data);
        pthread_mutex_unlock(&mutex);
        if(data[0] == 3.1415926 && data[1] == 3.1415926)
        {
            break;
        }
        else
        {
            send(sockfd,"Get!", 4, 0);
        }
    }

    close(sockfd);
    pthread_mutex_lock(&mutex);
    my_queue->len = -1;
    pthread_mutex_unlock(&mutex);
    printf("end TCP_client!!!\n");
}

void* MW_ASR(void* argv)
{
    queue* my_queue = (queue*)argv;
    double* data = NULL;
    double double[8]

    int sampling_rate = 500;
    int channels = 8;
    int cutoff = 5;
    double window_len = 0.5;
    double win_overlap = 0.66;
    double sig_window[channels][(int)window_len*sampling_rate];
    int sig_window_len = 0;
    double* filter_A_value = AB_filter(sampling_rate, 'A');
    double* filter_B_value = AB_filter(sampling_rate, 'B');
    ASR_PSW my_ASR = create_ASR(cutoff, window_len, win_overlap, sampling_rate, channels, filter_A_value, filter_B_value);
    my_ASR.data_length = window_len*sampling_rate;

    while(1)
    {
        if(!queue_is_empty(my_queue) && my_queue->len > 0)
        {
            pthread_mutex_lock(&mutex);
            data = pop_queue(my_queue);
            pthread_mutex_unlock(&mutex);
            if(data[0] == 3.1415926 && data[1] == 3.1415926) // session end signal
                break;
            if(sig_window_len >= window_len*sampling_rate)
            {

            }
            for(int i=0; i<channels; i++)
            {
                sig_window[i][sig_window_len] = data[i];
            }
            // printf("%d:\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n\n--------------------\n\n", i, data[0],data[1],data[2],data[3],data[4],data[5],data[6],data[7]);
        }
        else if(my_queue->len < 0)
            break;
    }

    queue_delete(my_queue);
    printf("end MW_ASR!!\n");
}

int cal_IE(double sig_window[][])
{

}

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <string.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

#define SOCKET_ERROR        -1
#define BUFFER_SIZE         100
#define MESSAGE             "This is the message I'm sending back and forth"
#define QUEUE_SIZE          5
#define HOST_NAME_SIZE      21

/* Common block to hold socket IDs  */
extern struct{
   int hsocket[8];
   int nhostport[8];
} skt_;


/*________________________________________________________________________*/

/* SOCKET SERVER OPEN*/
socketserveropen_()
{
    
    int hServerSocket;  /* handle to socket */
    struct hostent* pHostInfo;   /* holds info about a machine */
    struct sockaddr_in Address; /* Internet socket address stuct */
    int nAddressSize=sizeof(struct sockaddr_in);
    
    printf("\nMaking socket");
    /* make a socket */
     hServerSocket=socket(AF_INET,SOCK_STREAM,0);
   if(hServerSocket == SOCKET_ERROR)
    {
    
         printf("\nCould not make a socket.");
	return 0; 
    }
    
    
    /* fill address struct */
    Address.sin_addr.s_addr=INADDR_ANY;
    Address.sin_port=htons(skt_.nhostport[7]);
    Address.sin_family=AF_INET;

    printf("\nBinding to port %d",skt_.nhostport[7]);
    /* bind to a port */
    /*if(bind(hServerSocket,(struct sockaddr*)&Address,sizeof(Address)) 
                        == SOCKET_ERROR)*/
			
     while(bind(hServerSocket,(struct sockaddr*)&Address,sizeof(Address)) 
                        == SOCKET_ERROR){
	printf("\nCould not connect to host. Port \n", skt_.nhostport[7]);
           /*return 0;*/
	  skt_.nhostport[7]=skt_.nhostport[7]+1; 
	Address.sin_port=htons(skt_.nhostport[7]+1);
     }

    /*  get port number */
    getsockname( hServerSocket, (struct sockaddr *) &Address,(socklen_t *)&nAddressSize);
    printf("opened socket as fd (%d) on port (%d) for stream i/o\n",hServerSocket, 
    ntohs(Address.sin_port));

        printf("Server\n\
              sin_family        = %d\n\
              sin_addr.s_addr   = %d\n\
              sin_port          = %d\n"
              , Address.sin_family
              , Address.sin_addr.s_addr
              , ntohs(Address.sin_port)
            );
	    
            /* get port number to send it to main program */ 
            skt_.nhostport[7]=ntohs(Address.sin_port);


    printf("Making a listen queue of %d elements",QUEUE_SIZE);
    /* establish listen queue */
    if(listen(hServerSocket,QUEUE_SIZE) == SOCKET_ERROR)
    {
        printf("\nCould not listen\n");
        return 0;
    }

        printf("\nWaiting for a connection\n");

        skt_.hsocket[7]=accept(hServerSocket,(struct sockaddr*)&Address,(socklen_t *)&nAddressSize);

	printf("\nGot a connection");
 
        printf("\nSocketumber:\" %d \n",skt_.hsocket[7]);


}

	
/*________________________________________________________________________*/

/* SOCKET SERVER_CONNECT*/
socketserverconnect_()
{
    int hServerSocket;  /* handle to socket */
    struct hostent* pHostInfo;   /* holds info about a machine */
    struct sockaddr_in Address; /* Internet socket address stuct */
    int nAddressSize=sizeof(struct sockaddr_in);
 	
        /* get the connected socket */
        skt_.hsocket[7]=accept(hServerSocket,(struct sockaddr*)&Address,(socklen_t *)&nAddressSize);

	printf("\nGot a connection");
 
        printf("\nSocketumber:\" %d \n",skt_.hsocket[7]);
}





/*________________________________________________________________________*/

/* SOCKET CLIENT OPEN*/

socketclientopen_(hostportfc,strHostNamef)
int *hostportfc;
char strHostNamef[HOST_NAME_SIZE];
{
    struct hostent *pHostInfo;   /* holds info about a machine */
    struct sockaddr_in Address;  /* Internet socket address stuct */
    long nHostAddress;
    char pBuffer[BUFFER_SIZE];
    unsigned nReadAmount;
    int nhostportl;
    char strHostName[HOST_NAME_SIZE];
    
/* Copy input parameters to new variables */
    nhostportl=*hostportfc;
    strcpy(strHostName,strHostNamef);
    
    printf("\nMaking a socket");
    /* make a socket */
    skt_.hsocket[7]=socket(AF_INET,SOCK_STREAM,IPPROTO_TCP);
    if(skt_.hsocket[7] == SOCKET_ERROR)
    {
        printf("\nCould not make a socket\n");
        return 0;
    }


    /* get IP address from name */
    pHostInfo=gethostbyname(strHostName);
    /* copy address into long */
    memcpy(&nHostAddress,pHostInfo->h_addr,pHostInfo->h_length);

    /* fill address struct */
    Address.sin_addr.s_addr=nHostAddress;
    Address.sin_port=htons(nhostportl);
    Address.sin_family=AF_INET;

    printf("\nConnecting to %s on port %d:",strHostName,nhostportl);

    /* connect to host */
    if(connect(skt_.hsocket[7],(struct sockaddr*)&Address,sizeof(Address)) 
       == SOCKET_ERROR)
    {
        printf("\nCould not connect to host\n");
        return 0;
    }
}




/*________________________________________________________________________*/

/* SOCKET WRITE*/

socketwrite_(outvar, socketid)
float *outvar;
int *socketid;
{
char dummy1[100];
int stringlenght;

/*	usleep(6000000);  */

stringlenght=sprintf(dummy1,"%10.3e",*outvar);
stringlenght=stringlenght+1;
/*	printf("Sending [%s] with %d characters to Socket %d\n",dummy1,stringlenght, *socketid);*/
        /* write variable to socket.*/
        write(*socketid,dummy1,stringlenght);
}




/*________________________________________________________________________*/

/* SOCKET READ*/

socketread_(invar, socketid)
float *invar;
int *socketid;
{
char dummy1[BUFFER_SIZE];
char *pEnd;
double dummy2;
/*	usleep(6000000);    */

/*	printf("Reading sockets on c: %d.\n",*socketid); */

        /* read variable to socket.*/
        read(*socketid,dummy1,BUFFER_SIZE);

/*	printf("Variable received character format\"%s\n",dummy1); */
	dummy2=strtod(dummy1,&pEnd);
/*	printf("Converted to double precision \"%10.3e\n",dummy2); */
	*invar=dummy2;
/*	printf("Converted to floating variable \"%10.3f\n",*invar); */

}




/*________________________________________________________________________*/

/* SOCKET CLOSE*/

socketclose_(int socketid)
{
        close(socketid);
}

#include "stdint.h"
#include <stdlib.h>
#include "string.h"

struct bucket;

extern struct bucket bucket;

/*starts the first bucket*/
void *init(struct bucket, ssize_t capacity );


/* fills the bucket */
void *fill(struct bucket *bucket, void* data);

/*emptys the bucket*/
void *empty(struct bucket *bucket);

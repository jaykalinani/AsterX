This sets up a docker environment for developing CactusAMReX.

To use the docker environment, first run:

```bash
docker-compose build
```

This will create the local image for you to use.

Next, you need to start your image running. To do this, run:

```bash
docker-compose up -d
```

Once the image is running, you can get a shell inside it using this command:

```bash
docker exec -it carpetx bash
```

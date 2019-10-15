This sets up a docker environment for developing CactusAMReX.

To use the docker environment, first run:

  docker-compose build

This will create the local image for you to use.

Next, you need to start your image running. To do this, run:

  docker-compose up -d

Once the image is running, you can get a shell inside it using this command:

  docker-compose exec -it cacrex bash

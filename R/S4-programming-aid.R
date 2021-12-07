transfer_slot_content <- function(recipient, donor, skip = character(0), verbose = TRUE){

  snames_rec <- methods::slotNames(recipient)
  snames_don <- methods::slotNames(donor)

  for(snr in snames_rec){

    if(snr %in% snames_don & !snr %in% skip){

      give_feedback(
        msg = glue::glue("Transferring content of slot '{snr}'."),
        with.time = FALSE,
        verbose = verbose
      )

      recipient <-
        base::tryCatch({

          methods::slot(recipient, name = snr) <-
            methods::slot(donor, name = snr)

          recipient

        }, error = function(error){

          give_feedback(msg = error$message, verbose = verbose, with.time = FALSE)

          recipient

        })

    }

  }

  return(recipient)

}
